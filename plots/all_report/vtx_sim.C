// ==============================================================================================================
// vtx_sim.C -- copy (2026-11-14, integration of the event-vertex efficiency into report_all / report_exp) of
//   /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_verify/vv_sim.C
//   (study of the eta efficiency of the mandatory event-vertex cut vz > -500 mm; verifier ve_verify, own code).
// What it does: simulated eps_vtx per truth class C2-C5 (same base selection).
// Changes against the original: file / function / include names (vv_* -> vtx_*), and the text outputs go to
//   /hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/vtx_*.txt instead of the session scratchpad. Numbers, cuts and fits are unchanged.
// The original outputs are copied to docs/latex/tables_all/vtx_*.txt; the notes to docs/latex/notes_all/ve_*_note.md.
// Run (from /home/przygoda/HADES/FAT/FAT; vtx_fill.C first, then the analysis macros, at most 6 processes at once):
//   B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic", kTRUE)'
//   root -l -b -q -e "$B" 'plots/all_report/vtx_sim.C+'
// ==============================================================================================================
// vtx_sim.C -- ve_verify: eps_vtx(eta) in the simulation, own truth tagging (not sim_truth.h, not the ve_sim_cause macros)
// root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic",kTRUE)' 'vtx_sim.C+'
// Rows: simulated PAT files ppee_eta_all.root (s_00) and ppee_std_all.root (s_01), tree PipPimEpEm_ID, isBest >= 0,
// (trigbit & 8192) != 0. Truth (GEANT ids: e+ 2, e- 3, pi+ 8, pi- 9, pi0 7, eta 17, gamma 1; process 5 decay, 6 pair prod.):
//   C5 (eta set): the four tracks have the right ids, are direct products (parentid -1, process 0) of the forced eta decay
//       (geninfo == 1702030848.f, the float of code 1702030809) and share geninfo2 (same decay instance);
//   C4 (std set): pi+ and pi- with parent eta (17), process 5, same creation vertex; e+ and e- from pair production (6),
//       same creation vertex, grandparent eta (17); the lepton and pion generator codes (geninfo, geninfo1, geninfo2) equal;
//   C2 (std set): pions as C4; leptons with parent pi0 (7), process 5, grandparent eta, same vertex; codes equal;
//   C3 (std set): pions as C4; leptons pair production, same vertex, grandparent pi0 (7); codes equal.
#include <cstdio>
#include <cmath>
#include <vector>
#include "TChain.h"
#include "TString.h"

void vtx_sim() {
  FILE* o = fopen("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/vtx_sim.txt", "w");
  const float kEtaPPee = 1702030848.f;
  // counters [class][selection][pass]: selections 0 all M4 (cor_ok), 1 peak window 0.5333-0.5587, 2 bump window 0.4576-0.5255,
  //                                    3 data peak window x 1.0057 = 0.5367-0.5573 ;  pass: vz>-500 / vz>-160 / total
  double cnt[6][4][3] = {}; long long nocor[6][2] = {};
  const char* cname[6] = {"C5 eta->pi+pi-e+e- (eta set)", "C4 eta->pi+pi-gamma conv (std)", "C2 bump pi0 Dalitz (std)",
                          "C3 bump conversion (std)", "C4+C5", "C2+C3"};
  // lost rows: details of C5 rows without vertex
  long long lostNwin[5] = {};
  for (int set = 0; set < 2; ++set) {
    TChain ch("t"); ch.Add(Form("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/rows/s_%02d.root", set));
    Int_t cor, tid[4], tpar[4], tgp[4], tproc[4]; Float_t vz, m4, g0[4], g1[4], g2[4], tvx[4], tvy[4], tvz[4], z[4], r[4];
    ch.SetBranchAddress("cor", &cor); ch.SetBranchAddress("vz", &vz); ch.SetBranchAddress("m4", &m4);
    ch.SetBranchAddress("tid", tid); ch.SetBranchAddress("tpar", tpar); ch.SetBranchAddress("tgp", tgp); ch.SetBranchAddress("tproc", tproc);
    ch.SetBranchAddress("g0", g0); ch.SetBranchAddress("g1", g1); ch.SetBranchAddress("g2", g2);
    ch.SetBranchAddress("tvx", tvx); ch.SetBranchAddress("tvy", tvy); ch.SetBranchAddress("tvz", tvz);
    ch.SetBranchAddress("z", z); ch.SetBranchAddress("r", r);
    const Long64_t n = ch.GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
      ch.GetEntry(i);
      // slots: 0 e+, 1 e-, 2 pi+, 3 pi-
      const bool ids = tid[0] == 2 && tid[1] == 3 && tid[2] == 8 && tid[3] == 9;
      if (!ids) continue;
      auto sameV = [&](int a, int b) { return tvx[a] == tvx[b] && tvy[a] == tvy[b] && tvz[a] == tvz[b]; };
      auto sameC = [&](int a, int b) { return g0[a] == g0[b] && g1[a] == g1[b] && g2[a] == g2[b]; };
      int c = -1;
      if (set == 0) {
        bool ok = true;
        for (int k = 0; k < 4; ++k) if (!(g0[k] == kEtaPPee && tpar[k] == -1 && tproc[k] == 0 && g2[k] == g2[0])) ok = false;
        if (ok) c = 0;
      } else {
        const bool pions = tpar[2] == 17 && tpar[3] == 17 && tproc[2] == 5 && tproc[3] == 5 && sameV(2, 3) && sameC(2, 3);
        if (pions && sameC(0, 2) && sameC(1, 2)) {
          if (tproc[0] == 6 && tproc[1] == 6 && sameV(0, 1) && tgp[0] == 17 && tgp[1] == 17) c = 1;
          else if (tpar[0] == 7 && tpar[1] == 7 && tproc[0] == 5 && tproc[1] == 5 && sameV(0, 1) && tgp[0] == 17 && tgp[1] == 17) c = 2;
          else if (tproc[0] == 6 && tproc[1] == 6 && sameV(0, 1) && tgp[0] == 7 && tgp[1] == 7) c = 3;
        }
      }
      if (c < 0) continue;
      const int cc[2] = {c, c <= 1 ? 4 : 5};
      for (int q = 0; q < 2; ++q) {
        const int k = cc[q];
        if (!cor) { nocor[k][0]++; if (vz > -500) nocor[k][1]++; continue; }
        bool w[4] = {true, m4 >= 0.5333 && m4 < 0.5587, m4 >= 0.4576 && m4 < 0.5255, m4 >= 0.53671 && m4 < 0.55732};
        for (int s = 0; s < 4; ++s) if (w[s]) { cnt[k][s][2] += 1; if (vz > -500) cnt[k][s][0] += 1; if (vz > -160) cnt[k][s][1] += 1; }
      }
      if (c == 0 && cor && !(vz > -500)) {
        int nw = 0; for (int k = 0; k < 4; ++k) if (z[k] > -170 && z[k] < -30 && std::fabs(r[k]) < 20) ++nw; lostNwin[std::min(nw, 4)]++;
      }
    }
  }
  auto P = [&](const char* s) { printf("%s", s); fprintf(o, "%s", s); };
  P("# vtx_sim: eps_vtx = rows with vz > -500 / all rows (isBest >= 0, PT3 bit, cor_ok), own truth tagging, unweighted\n");
  P("# windows: sim-centred peak 0.5333-0.5587, bump 0.4576-0.5255, data peak window x 1.0057 = 0.5367-0.5573\n");
  const char* wn[4] = {"all M4", "peak(sim)", "bump(sim)", "peak(data x1.0057)"};
  for (int k = 0; k < 6; ++k) {
    P(Form("%s\n", cname[k]));
    for (int s = 0; s < 4; ++s) {
      const double N = cnt[k][s][2]; if (N <= 0) continue;
      const double e = cnt[k][s][0] / N, ee = std::sqrt(std::max(e * (1 - e), 1.0 / N) / N);
      const double e2 = cnt[k][s][1] / N;
      // one-sided 95% lower limit when all pass: 0.05^(1/N)
      P(Form("   %-20s rows %8.0f  eps(vz>-500) %.5f +- %.5f  (lost %.0f; 95%%CL lower %.4f)  eps(vz>-160) %.5f\n", wn[s], N, e, ee,
             N - cnt[k][s][0], cnt[k][s][0] == N ? std::pow(0.05, 1.0 / N) : e - 1.645 * ee, e2));
    }
    P(Form("   rows with cor_ok == 0: %lld (vz>-500: %lld)\n", nocor[k][0], nocor[k][1]));
  }
  P(Form("C5 rows without vertex, own tracks in the finder window (0..4): %lld %lld %lld %lld %lld\n", lostNwin[0], lostNwin[1], lostNwin[2],
         lostNwin[3], lostNwin[4]));
  fclose(o);
}
