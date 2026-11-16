// ==============================================================================================================
// vtx_elastic.C -- copy (2026-11-14, integration of the event-vertex efficiency into report_all / report_exp) of
//   /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_sim_cause/elastic.C
//   (study of the eta efficiency of the mandatory event-vertex cut vz > -500 mm; package ve_sim_cause).
// What it does: pp-elastic content of 2-track events with / without a vertex (is the upstream population hydrogen?).
// Changes against the original: only the file and function names (-> vtx_*); the code is unchanged (it writes the ROOT file
//   given as argument; the summed result of the original run is ve_sim_cause/elsum_all.txt, quoted in notes_all/ve_sim_cause_note.md).
// The original outputs are copied to docs/latex/tables_all/vtx_*.txt; the notes to docs/latex/notes_all/ve_*_note.md.
// Run (from /home/przygoda/HADES/FAT/FAT; vtx_fill.C first, then the analysis macros, at most 6 processes at once):
//   B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic", kTRUE)'
//   for d in 036 049 062; do root -l -b -q "plots/all_report/vtx_elastic.C+(\"/hdd2/przygoda/HADES/pp45/exp/had_${d}_p1.root\",-1,\"/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_sim_cause/elastic_full_${d}_p1.root\")"; done
// ==============================================================================================================
// elastic.C -- is the upstream population (track z < -170 mm, no event vertex) hydrogen?
// pp elastic test on the hadron PAT tree PP_ID (two tracks under the proton hypothesis):
// events with exactly two tracks (totalmult == 2), isBest == 1; coplanarity |dphi - 180| and
// tan(theta1) tan(theta2) = 1/gamma_cm^2 = 0.2943 at T = 4.5 GeV, in bins of zbar = (z1 + z2)/2
// (track z = closest approach to the beam axis), with and without event vertex.
// Usage: root -l -b -q 'vtx_elastic.C+("had.root", maxRows, "out.root")'
#include <cmath>
#include <cstdio>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"

void vtx_elastic(const char* fname, long long maxRows, const char* out) {
  TFile* f = TFile::Open(fname);
  TTree* t = (TTree*)f->Get("PP_ID");
  t->SetBranchStatus("*", 0);
  float isBest, vz, mult, th1, th2, ph1, ph2, z1, z2, r1, r2, p1, p2, st, tb;
  auto B = [&](const char* n, float* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
  B("isBest", &isBest); B("eVertReco_z", &vz); B("totalmult", &mult); B("evtStart", &st); B("trigbit", &tb);
  B("p1_theta", &th1); B("p2_theta", &th2); B("p1_phi", &ph1); B("p2_phi", &ph2);
  B("p1_z", &z1); B("p2_z", &z2); B("p1_r", &r1); B("p2_r", &r2); B("p1_p", &p1); B("p2_p", &p2);
  const int nz = 7;
  const double ze[nz + 1] = {-700, -300, -240, -200, -170, -160, -90, 0};
  TFile fo(out, "RECREATE");
  TH1D* htt[2][nz]; TH1D* hdp[2][nz];
  for (int v = 0; v < 2; ++v)
    for (int i = 0; i < nz; ++i) {
      htt[v][i] = new TH1D(Form("tt_%s_z%d", v ? "nov" : "vtx", i), Form("tan th1 tan th2, |dphi-180|<3, zbar [%g,%g)", ze[i], ze[i + 1]), 100, 0, 1);
      hdp[v][i] = new TH1D(Form("dp_%s_z%d", v ? "nov" : "vtx", i), Form("dphi, |tt-0.294|<0.04, zbar [%g,%g)", ze[i], ze[i + 1]), 90, 135, 225);
    }
  TH1D hz0("zbar_nov", "zbar, 2-track events without vertex", 140, -700, 0), hz1("zbar_vtx", "zbar, with vertex", 140, -700, 0);
  TH1D hze0("zbar_el_nov", "zbar, elastic candidates, no vertex", 140, -700, 0), hze1("zbar_el_vtx", "zbar, elastic candidates, vertex", 140, -700, 0);
  long long n = t->GetEntries(); if (maxRows > 0 && maxRows < n) n = maxRows;
  long long nsel = 0, nnov = 0;
  for (long long e = 0; e < n; ++e) {
    t->GetEntry(e);
    if (isBest < 0.5 || mult != 2) continue;
    ++nsel;
    const bool nov = vz <= -999.f; if (nov) ++nnov;
    const double zb = 0.5 * (z1 + z2);
    (nov ? hz0 : hz1).Fill(zb);
    double dp = std::fabs(ph1 - ph2); if (dp > 360) dp -= 360;
    const double tt = std::tan(th1 * TMath::DegToRad()) * std::tan(th2 * TMath::DegToRad());
    int iz = -1; for (int i = 0; i < nz; ++i) if (zb >= ze[i] && zb < ze[i + 1]) iz = i;
    if (iz < 0) continue;
    if (std::fabs(dp - 180) < 3) htt[nov][iz]->Fill(tt);
    if (std::fabs(tt - 0.2943) < 0.04) hdp[nov][iz]->Fill(dp);
    if (std::fabs(dp - 180) < 3 && std::fabs(tt - 0.2943) < 0.04) (nov ? hze0 : hze1).Fill(zb);
  }
  printf("file %s rows read %lld, 2-track isBest events %lld, without vertex %lld (%.4f)\n", fname, n, nsel, nnov, (double)nnov / nsel);
  // peak / side-band summary: elastic peak |tt - 0.2943| < 0.04 vs side bands 0.15-0.23 and 0.36-0.44 (same width x2)
  for (int v = 0; v < 2; ++v)
    for (int i = 0; i < nz; ++i) {
      TH1D* h = htt[v][i];
      const double pk = h->Integral(h->FindBin(0.2543 + 1e-6), h->FindBin(0.3343 - 1e-6));
      const double sb = 0.5 * (h->Integral(h->FindBin(0.15 + 1e-6), h->FindBin(0.23 - 1e-6)) + h->Integral(h->FindBin(0.36 + 1e-6), h->FindBin(0.44 - 1e-6)));
      TH1D* d = hdp[v][i];
      const double dpk = d->Integral(d->FindBin(177.01), d->FindBin(182.99));
      const double dsb = 0.5 * (d->Integral(d->FindBin(160.01), d->FindBin(165.99)) + d->Integral(d->FindBin(194.01), d->FindBin(199.99)));
      printf("  %s zbar [%5.0f,%5.0f): coplanar rows %8.0f | tt peak %8.0f side %8.0f -> excess %8.0f (peak/side %.2f) | dphi peak %8.0f side %8.0f (%.2f)\n",
             v ? "no vertex" : "vertex   ", ze[i], ze[i + 1], h->GetEntries(), pk, sb, pk - sb, sb > 0 ? pk / sb : 0, dpk, dsb, dsb > 0 ? dpk / dsb : 0);
    }
  fo.Write(); fo.Close();
}
