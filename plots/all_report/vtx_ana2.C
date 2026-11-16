// ==============================================================================================================
// vtx_ana2.C -- copy (2026-11-14, integration of the event-vertex efficiency into report_all / report_exp) of
//   /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_verify/vv_ana2.C
//   (study of the eta efficiency of the mandatory event-vertex cut vz > -500 mm; verifier ve_verify, own code).
// What it does: N-only fits, the Gaussian content of the CB of Nv (e+e+ dip), pi0 tag in Nv, omega-region test.
// Changes against the original: file / function / include names (vv_* -> vtx_*), and the text outputs go to
//   /hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/vtx_*.txt instead of the session scratchpad. Numbers, cuts and fits are unchanged.
// The original outputs are copied to docs/latex/tables_all/vtx_*.txt; the notes to docs/latex/notes_all/ve_*_note.md.
// Run (from /home/przygoda/HADES/FAT/FAT; vtx_fill.C first, then the analysis macros, at most 6 processes at once):
//   B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic", kTRUE)'
//   root -l -b -q -e "$B" 'plots/all_report/vtx_ana2.C+'
// ==============================================================================================================
// vtx_ana2.C -- ve_verify: further tests of the eta content of the no-vertex rows (own code)
//   A. where does S(Nv) come from: fit of N alone vs N - CB; the CB's own content at m_eta; like-sign rows per period
//   B. pi0 peak in M(e+e-gamma) in V and Nv (is the photon tag functional without an event vertex?)
//   C. the structure near M4 = 0.73-0.75 (omega region) in V and Nv: another physics-content test
//   D. T split by the pion crossing class; T and E per period
// root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic",kTRUE)' 'vtx_ana2.C+'
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "TChain.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TString.h"
#include "vtx_lib.h"
using namespace vvl;
#include "vtx_load.h"

static FILE* gO = nullptr;
#define P(...) do { printf(__VA_ARGS__); if (gO) fprintf(gO, __VA_ARGS__); } while (0)
static int Per(int d) { return d <= 46 ? 0 : d <= 48 ? 1 : d <= 53 ? 2 : 3; }

void vtx_ana2() {
  gO = fopen("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/vtx_ana2.txt", "w");
  std::vector<R> d = Load("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/rows/d_*.root");
  const char* pn[4] = {"I", "IIA", "IIB", "III"};
  std::vector<std::string> nm = {"V", "Nv"};
  std::vector<Sel> s = {[](const R& r) { return r.vz > -500; }, [](const R& r) { return r.vz <= -500; }};
  for (int p = 0; p < 4; ++p) {
    nm.push_back(Form("V,per%s", pn[p])); s.push_back([p](const R& r) { return r.vz > -500 && Per(r.day) == p; });
    nm.push_back(Form("Nv,per%s", pn[p])); s.push_back([p](const R& r) { return r.vz <= -500 && Per(r.day) == p; });
  }
  std::vector<Spec> sp = Fill(d, s);

  P("# A. S from N - CB (P_main), from N alone (no CB; pol2 and pol3), and the CB's own Gaussian content at m_eta\n");
  P("  %-12s %8s %6s | %8s %6s %8s %6s | %8s %6s | %8s %6s %8s %6s\n", "class", "S(N-CB)", "e", "S(N)p2", "e", "S(N)p3", "e", "S_CB", "e", "S(N++)", "e", "S(N--)", "e");
  for (size_t i = 0; i < s.size(); ++i) {
    PCfg c0; Fit a = Peak(sp[i], c0);
    PCfg c5; c5.what = 5; Fit b = Peak(sp[i], c5);
    PCfg c5b; c5b.what = 5; c5b.ord = 3; Fit b3 = Peak(sp[i], c5b);
    PCfg c6; c6.what = 6; Fit cb = Peak(sp[i], c6);
    PCfg c1; c1.what = 1; Fit pp = Peak(sp[i], c1);
    PCfg c2; c2.what = 2; Fit mm = Peak(sp[i], c2);
    P("  %-12s %8.0f %6.0f | %8.0f %6.0f %8.0f %6.0f | %8.0f %6.0f | %8.0f %6.0f %8.0f %6.0f\n", nm[i].c_str(), a.S(), a.eS(), b.S(), b.eS(), b3.S(),
      b3.eS(), cb.S(), cb.eS(), pp.S(), pp.eS(), mm.S(), mm.eS());
  }
  // per-period ratio of Nv eta per Nv us row
  P("  eta per 1e3 us rows (S(N-CB) / rows): ");
  for (size_t i = 0; i < s.size(); ++i) { PCfg c0; Fit a = Peak(sp[i], c0); P(" %s %.3f+-%.3f;", nm[i].c_str(), 1e3 * a.S() / sp[i].nsel[0], 1e3 * a.eS() / sp[i].nsel[0]); }
  P("\n");

  // ---- B. pi0 in M(ee gamma) ----
  P("\n# B. pi0 peak in M(e+e-gamma) (rows with a good photon, M(ee) < 0.10; N - CB, 2.5 MeV bins, Gaussian + pol2 in 0.08-0.20)\n");
  {
    TH1D* h[2][3];
    for (int c = 0; c < 2; ++c) for (int k = 0; k < 3; ++k) h[c][k] = new TH1D(Form("meeg_%d_%d", c, k), "", 120, 0.0, 0.30);
    for (const R& r : d) {
      if (!r.cor || r.ng <= 0 || !(r.mee < 0.10) || r.meeg < 0) continue;
      h[r.vz > -500 ? 0 : 1][r.smp]->Fill(r.meeg);
    }
    for (int c = 0; c < 2; ++c) {
      TH1D* sig = (TH1D*)h[c][0]->Clone(Form("sig%d", c));
      for (int b = 1; b <= 120; ++b) {
        const double a = h[c][1]->GetBinContent(b), m = h[c][2]->GetBinContent(b), n = h[c][0]->GetBinContent(b);
        sig->SetBinContent(b, n - 2 * std::sqrt(a * m)); sig->SetBinError(b, std::sqrt(std::max(n + a + m, 1.0)));
      }
      TF1 f(Form("f%d", c), "[0]*0.0025/(sqrt(2*pi)*[2])*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x+[5]*x*x", 0.08, 0.20);
      f.SetParameters(sig->Integral(sig->FindBin(0.12), sig->FindBin(0.15)), 0.135, 0.010, 0, 0, 0);
      sig->Fit(&f, "QR0");
      P("  %-3s pi0: yield %.0f +- %.0f, mean %.4f +- %.4f, sigma %.4f +- %.4f GeV/c2, chi2/ndf %.1f/%d; N - CB in 0.119-0.173: %.0f of N %.0f\n",
        c == 0 ? "V" : "Nv", f.GetParameter(0), f.GetParError(0), f.GetParameter(1), f.GetParError(1), f.GetParameter(2), f.GetParError(2),
        f.GetChisquare(), f.GetNDF(), sig->Integral(sig->FindBin(0.1195), sig->FindBin(0.1725)), h[c][0]->Integral(h[c][0]->FindBin(0.1195), h[c][0]->FindBin(0.1725)));
    }
  }

  // ---- C. structure near 0.73-0.75 ----
  P("\n# C. M4 structure in the omega region: frozen Gaussian (5.4 and 8 MeV) + pol2 at the positions of the largest V pulls\n");
  for (double mu : {0.7325, 0.7400}) for (double sg : {0.0054, 0.0080}) {
    PCfg c; c.mu = mu; c.sg = sg; c.lo = mu - 0.0214; c.hi = mu + 0.0661;
    Fit a = Peak(sp[0], c), b = Peak(sp[1], c); double e; const double v = EpsV(a.S(), a.eS(), b.S(), b.eS(), e);
    P("  mu %.4f sigma %.1f MeV: S_V %.0f +- %.0f, S_Nv %.0f +- %.0f, eps %.3f +- %.3f, ratio Nv/V %.3f +- %.3f\n", mu, sg * 1e3, a.S(), a.eS(),
      b.S(), b.eS(), v, e, b.S() / a.S(), std::sqrt(std::pow(b.eS() / a.S(), 2) + std::pow(b.S() * a.eS() / (a.S() * a.S()), 2)));
  }
  // broad window content of N - CB in 0.70-0.76 relative to a pol2 fitted in 0.62-0.70 + 0.80-0.90 (side bands)
  {
    for (int i = 0; i < 2; ++i) {
      // linear GLS on side bands only
      std::vector<int> idx; for (int b = 0; b < kNb; ++b) { const double x = sp[i].x(b); if ((x >= 0.62 && x < 0.70) || (x >= 0.80 && x < 0.90)) idx.push_back(b); }
      std::vector<std::vector<double>> F; std::vector<double> y, add;
      for (int b : idx) { const double t = (sp[i].x(b) - 0.76) / 0.14; F.push_back({1.0, t, t * t, t * t * t}); y.push_back(sp[i].n[0][b] - sp[i].cb[b]); add.push_back(sp[i].cb[b]); }
      // GLS needs S at index 0; use the constant term as "S" index (only the convergence test uses it)
      Fit f = GLS(sp[i], idx, F, y, add, true);
      double ex = 0, var = 0;
      for (int b = 0; b < kNb; ++b) { const double x = sp[i].x(b); if (x < 0.70 || x >= 0.76) continue; const double t = (x - 0.76) / 0.14;
        ex += sp[i].n[0][b] - sp[i].cb[b] - (f.p[0] + f.p[1] * t + f.p[2] * t * t + f.p[3] * t * t * t); var += sp[i].n[0][b] + sp[i].cb[b] * 0.0; }
      P("  %s: excess of N - CB over a pol3 side-band fit (0.62-0.70, 0.80-0.90) in 0.70-0.76: %.0f +- %.0f (stat of N only)\n", i == 0 ? "V " : "Nv", ex, std::sqrt(var));
    }
  }
  fclose(gO);
}
