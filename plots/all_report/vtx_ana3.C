// ==============================================================================================================
// vtx_ana3.C -- copy (2026-11-14, integration of the event-vertex efficiency into report_all / report_exp) of
//   /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_verify/vv_ana3.C
//   (study of the eta efficiency of the mandatory event-vertex cut vz > -500 mm; verifier ve_verify, own code).
// What it does: S(V), S(Nv) with a smooth like-sign CB (12 variants), CB in 20 MeV bins, wide Gaussians.
// Changes against the original: file / function / include names (vv_* -> vtx_*), and the text outputs go to
//   /hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/vtx_*.txt instead of the session scratchpad. Numbers, cuts and fits are unchanged.
// The original outputs are copied to docs/latex/tables_all/vtx_*.txt; the notes to docs/latex/notes_all/ve_*_note.md.
// Run (from /home/przygoda/HADES/FAT/FAT; vtx_fill.C first, then the analysis macros, at most 6 processes at once):
//   B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic", kTRUE)'
//   root -l -b -q -e "$B" 'plots/all_report/vtx_ana3.C+'
// ==============================================================================================================
// vtx_ana3.C -- ve_verify: S(V), S(Nv) with a SMOOTH like-sign CB (no statistical structure of the like-sign rows at m_eta),
// with 20 MeV CB bins, and N-only fits with wider Gaussians (search for a smeared component)
// root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic",kTRUE)' 'vtx_ana3.C+'
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include "TChain.h"
#include "TString.h"
#include "vtx_lib.h"
using namespace vvl;
#include "vtx_load.h"

static FILE* gO = nullptr;
#define P(...) do { printf(__VA_ARGS__); if (gO) fprintf(gO, __VA_ARGS__); } while (0)

// smooth fit of a like-sign fine-bin spectrum with a polynomial of order ord in [lo, hi) (Poisson variance, iterated)
std::vector<double> Smooth(const Spec& sp, int k, double lo, double hi, int ord, double& chi2, int& ndf) {
  std::vector<int> idx; for (int i = 0; i < kNb; ++i) if (sp.x(i) >= lo && sp.x(i) < hi) idx.push_back(i);
  const double x0 = 0.5 * (lo + hi), dx = 0.5 * (hi - lo);
  std::vector<std::vector<double>> F; std::vector<double> y, add;
  for (int i : idx) { std::vector<double> row; const double t = (sp.x(i) - x0) / dx; double tp = 1;
    for (int j = 0; j <= ord; ++j) { row.push_back(tp); tp *= t; } F.push_back(row); y.push_back(sp.n[k][i]); add.push_back(0); }
  Fit f = GLS(sp, idx, F, y, add, false);
  chi2 = f.chi2; ndf = f.ndf;
  std::vector<double> out(kNb, 0.);
  for (int i = 0; i < kNb; ++i) { const double t = (sp.x(i) - x0) / dx; double tp = 1, v = 0; for (int j = 0; j <= ord; ++j) { v += f.p[j] * tp; tp *= t; } out[i] = v; }
  return out;
}

void vtx_ana3() {
  gO = fopen("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/vtx_ana3.txt", "w");
  std::vector<R> d = Load("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/rows/d_*.root");
  std::vector<std::string> nm = {"V", "Nv", "ALL"};
  std::vector<Sel> s = {[](const R& r) { return r.vz > -500; }, [](const R& r) { return r.vz <= -500; }, [](const R&) { return true; }};
  std::vector<Spec> sp = Fill(d, s);
  P("# A. smooth like-sign CB: N++ and N-- fitted with pol(ord) in [lo, hi); CB = 2 sqrt(f++ f--) per fine bin; P_main on N - CB\n");
  P("#    (no CB covariance: the smooth CB carries no bin-to-bin noise; its overall error is absorbed by pol2)\n");
  double SV[12], SN[12], eV[12], eN[12]; int nv = 0;
  for (int ord : {3, 4, 5}) for (double lo : {0.45, 0.48}) for (double hi : {0.66, 0.70}) {
    double res[3][2];
    for (int c = 0; c < 2; ++c) {
      double c1, c2; int n1, n2;
      std::vector<double> fp = Smooth(sp[c], 1, lo, hi, ord, c1, n1), fm = Smooth(sp[c], 2, lo, hi, ord, c2, n2);
      Spec q = sp[c];
      for (int i = 0; i < kNb; ++i) q.cb[i] = (fp[i] > 0 && fm[i] > 0) ? 2 * std::sqrt(fp[i] * fm[i]) : 0;
      for (int j = 0; j < q.nw; ++j) q.cbWv[j] = 0;   // no CB covariance
      PCfg pc; Fit f = Peak(q, pc);
      res[c][0] = f.S(); res[c][1] = f.eS();
      if (c == 1 && ord == 4 && lo == 0.45 && hi == 0.70) P("    (Nv, pol4 0.45-0.70: like-sign smooth fits chi2 %.1f/%d (++), %.1f/%d (--))\n", c1, n1, c2, n2);
    }
    double e; const double v = EpsV(res[0][0], res[0][1], res[1][0], res[1][1], e);
    P("  ord %d range %.2f-%.2f: S_V %6.0f +- %4.0f   S_Nv %6.0f +- %4.0f   eps %.3f +- %.3f\n", ord, lo, hi, res[0][0], res[0][1], res[1][0], res[1][1], v, e);
    SV[nv] = res[0][0]; SN[nv] = res[1][0]; eV[nv] = res[0][1]; eN[nv] = res[1][1]; ++nv;
  }
  double mV = 0, mN = 0; for (int i = 0; i < nv; ++i) { mV += SV[i]; mN += SN[i]; } mV /= nv; mN /= nv;
  double rN = 0; for (int i = 0; i < nv; ++i) rN += (SN[i] - mN) * (SN[i] - mN); rN = std::sqrt(rN / nv);
  P("  mean over %d smooth-CB variants: S_V %.0f, S_Nv %.0f (rms %.0f)\n", nv, mV, mN, rN);

  P("\n# B. CB in 20 MeV bins (8 fine bins), linear interpolation, P_main\n");
  {
    // rebuild with reb = 8 by summing pairs of 10 MeV CB bins (own implementation)
    for (int c = 0; c < 2; ++c) {
      Spec q = sp[c];
      const int nw2 = kNb / 8;
      std::vector<double> cb2(nw2), v2(nw2);
      for (int j = 0; j < nw2; ++j) { double a = 0, b = 0; for (int i = 8 * j; i < 8 * j + 8; ++i) { a += q.n[1][i]; b += q.n[2][i]; }
        cb2[j] = 2 * std::sqrt(a * b); v2[j] = a + b; }
      for (int i = 0; i < kNb; ++i) {
        const double t = (q.x(i) - (kLo + 0.5 * 8 * kBw)) / (8 * kBw); int j0, j1; double f;
        if (t <= 0) { j0 = j1 = 0; f = 0; } else if (t >= nw2 - 1) { j0 = j1 = nw2 - 1; f = 0; } else { j0 = (int)t; j1 = j0 + 1; f = t - j0; }
        q.cb[i] = ((1 - f) * cb2[j0] + f * cb2[j1]) / 8;
      }
      for (int j = 0; j < q.nw; ++j) q.cbWv[j] *= 0.5;   // approximate: CB variance per fine bin halves with twice wider bins
      PCfg pc; Fit f = Peak(q, pc);
      P("  %s: S = %.0f +- %.0f\n", nm[c].c_str(), f.S(), f.eS());
    }
  }

  P("\n# C. N-only fits (no CB; pol2 / pol3 in 0.5225-0.61) with wider Gaussians at 0.5439: a smeared eta component?\n");
  for (double sg : {0.0054, 0.0080, 0.0120}) for (int ord : {2, 3}) {
    PCfg c; c.what = 5; c.sg = sg; c.ord = ord;
    Fit a = Peak(sp[0], c), b = Peak(sp[1], c);
    PCfg c0; c0.sg = sg; c0.ord = ord; Fit a0 = Peak(sp[0], c0), b0 = Peak(sp[1], c0);
    P("  sigma %4.1f MeV pol%d: N only: S_V %6.0f +- %4.0f  S_Nv %6.0f +- %4.0f | N - CB: S_V %6.0f +- %4.0f  S_Nv %6.0f +- %4.0f\n", sg * 1e3, ord,
      a.S(), a.eS(), b.S(), b.eS(), a0.S(), a0.eS(), b0.S(), b0.eS());
  }
  fclose(gO);
}
