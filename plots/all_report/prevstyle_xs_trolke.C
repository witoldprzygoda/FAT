// ============================================================================
// prevstyle_xs_trolke.C  --  package prevstyle, section 10: TRolke (W. A. Rolke, A. M. Lopez, J. Conrad, NIM A 551 (2005) 493)
//                            90 % CL intervals of the X17 counts in one counting window, earlier-analysis style
//                            DEMONSTRATION (not a measurement; no result of the report uses it)
// ============================================================================
// Written 2026-11-16. New file. Called by prevstyle_results_xs.py (it writes the input list and reads the output).
// Input : text file, one line per configuration:  <id> <x = observed counts, integer> <B = expected background> <sigma_B>
//         (B and sigma_B from the side-band estimate: CB + k x side bands + eta-peak remainder, see prevstyle_results_xs.py)
// Output: stdout, one line per configuration and bounding option:  <id> <bounded> <lower> <upper>
// Model : model 5 of TRolke, SetGaussBkgKnownEff(x, B, sigma_B, 1.0): Poisson count, Gaussian background, known efficiency (1:
//         the interval is on the counts; the efficiency enters later in the BR conversion). Same call as trolke_x17.C.
// Run   : root -l -b -q 'plots/all_report/prevstyle_xs_trolke.C("<input file>")'
// ============================================================================
#include "TRolke.h"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

void prevstyle_xs_trolke(const char* infile) {
    std::ifstream in(infile);
    if (!in) { printf("ERROR cannot open %s\n", infile); return; }
    std::string line;
    int n = 0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream is(line);
        std::string id;
        long x = 0;
        double B = 0, sB = 0;
        if (!(is >> id >> x >> B >> sB)) continue;
        for (int bounded = 1; bounded >= 0; --bounded) {
            TRolke r(0.90);
            r.SetBounding(bounded);
            r.SetGaussBkgKnownEff((int)x, B, sB, 1.0);
            double lo = 0, hi = 0;
            r.GetLimits(lo, hi);
            printf("%s %d %.4f %.4f\n", id.c_str(), bounded, lo, hi);
        }
        ++n;
    }
    printf("# done %d configurations\n", n);
}
