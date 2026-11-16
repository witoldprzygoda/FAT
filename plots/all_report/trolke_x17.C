// trolke_x17.C -- TRolke (Rolke, Lopez, Conrad, NIM A 551 (2005) 493) on the frozen X17c counting window, for comparison
// with the 2D CLs limit. Run: root -l -b -q plots/all_report/trolke_x17.C > docs/latex/tables_all/trolke_x17.txt
// Inputs from x17d_results.md sect. 2.7: N = 6424 (exact count) or 6466 (overlap rule), B_pred = 6311.7 +- 1195.5 (frozen),
// BR = excess/(a_eta s1), a_eta = 0.3823, s1 = 2.0954e7.
#include "TRolke.h"
#include <cstdio>
void trolke_x17() {
  const double B = 6311.7, sB = 1195.5;   // frozen background prediction and its error
  const double scale = 1.0 / (0.3823 * 2.0954e7); // BR per row: BR = excess/(a_eta*s1)
  const int counts[2] = {6424, 6466};          // exact window count (rounded edges 6424); report value B+154
  for (int ic = 0; ic < 2; ++ic) {
    const int x = counts[ic];
    for (int bounded = 0; bounded < 2; ++bounded) {
      TRolke r(0.90);
      r.SetBounding(bounded);
      r.SetGaussBkgKnownEff(x, B, sB, 1.0);   // model 5: Gaussian background, known efficiency
      double lo = 0, hi = 0; r.GetLimits(lo, hi);
      printf("x=%d bounded=%d: N_X in [%.0f, %.0f] rows -> BR x BR in [%.2e, %.2e]\n", x, bounded, lo, hi, lo*scale, hi*scale);
    }
  }
}
