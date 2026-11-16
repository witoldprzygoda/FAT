// ==============================================================================================================
// vtx_k0s.C -- copy (2026-11-14, integration of the event-vertex efficiency into report_all / report_exp) of
//   /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_verify/vv_k0s.C
//   (study of the eta efficiency of the mandatory event-vertex cut vz > -500 mm; verifier ve_verify, own code).
// What it does: K0S -> pi+pi- in V and Nv (side check).
// Changes against the original: file / function / include names (vv_* -> vtx_*), and the text outputs go to
//   /hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/vtx_*.txt instead of the session scratchpad. Numbers, cuts and fits are unchanged.
// The original outputs are copied to docs/latex/tables_all/vtx_*.txt; the notes to docs/latex/notes_all/ve_*_note.md.
// Run (from /home/przygoda/HADES/FAT/FAT; vtx_fill.C first, then the analysis macros, at most 6 processes at once):
//   B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic", kTRUE)'
//   root -l -b -q -e "$B" 'plots/all_report/vtx_k0s.C+'
// ==============================================================================================================
// vtx_k0s.C -- ve_verify: quick K0S -> pi+pi- check in V and Nv (M(pi+pi-) of rows with M4 < 1.10; N - CB, CB per 10 MeV of M(pipi))
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include "TChain.h"
#include "TString.h"
#include "vtx_lib.h"
using namespace vvl;
#include "vtx_load.h"
void vtx_k0s() {
  FILE* o = fopen("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/vtx_k0s.txt", "w");
  std::vector<R> d = Load("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/rows/d_*.root");
  std::vector<Sel> s = {[](const R& r) { return r.vz > -500 && r.m4 < 1.10; }, [](const R& r) { return r.vz <= -500 && r.m4 < 1.10; }};
  std::vector<Spec> sp = Fill(d, s, [](const R& r) { return r.mpp; }, 0.485, 0.505);
  const char* nm[2] = {"V", "Nv"};
  for (double lo : {0.455, 0.46}) for (double hi : {0.545, 0.54}) for (double sg : {0.0060, 0.0064, 0.0075, 0.0089}) for (int ord : {2, 3}) {
    PCfg c; c.mu = 0.4950; c.sg = sg; c.lo = lo; c.hi = hi; c.ord = ord;
    Fit a = Peak(sp[0], c), b = Peak(sp[1], c); double e; const double v = EpsV(a.S(), a.eS(), b.S(), b.eS(), e);
    fprintf(o, "range %.3f-%.3f sigma %.1f MeV pol%d: Y_V %6.0f +- %4.0f (chi2 %.1f/%d)  Y_Nv %6.0f +- %4.0f (chi2 %.1f/%d)  eps %.3f +- %.3f\n", lo, hi, sg * 1e3, ord,
            a.S(), a.eS(), a.chi2, a.ndf, b.S(), b.eS(), b.chi2, b.ndf, v, e);
  }
  // Nv with its own wider width, V with nominal
  for (double sgN : {0.0060, 0.0075, 0.0089, 0.0100}) {
    PCfg cv; cv.mu = 0.4950; cv.sg = 0.0062; cv.lo = 0.455; cv.hi = 0.545;
    PCfg cn = cv; cn.sg = sgN; Fit a = Peak(sp[0], cv), b = Peak(sp[1], cn); double e; const double v = EpsV(a.S(), a.eS(), b.S(), b.eS(), e);
    fprintf(o, "V sigma 6.2, Nv sigma %.1f: Y_V %.0f +- %.0f, Y_Nv %.0f +- %.0f, eps %.3f +- %.3f\n", sgN * 1e3, a.S(), a.eS(), b.S(), b.eS(), v, e);
  }
  fclose(o);
}
