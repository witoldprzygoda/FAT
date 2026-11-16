// ============================================================================
// wp2_runs.C  --  WP2: per-run counts of S0 rows with / without a sector-2 track (time span of the deficit)
// ============================================================================
// Writes /hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_runs.csv:
//   run, day, then for us / pp / mm: n and variance of tot, s2, no2, A12, A34, B012, B3 (S0, all M4; variance = event-group)
//   A = no sector-2 track, B = >= 1 sector-2 track; 12 / 012 = fewer than 3 distinct non-2 sectors, 34 / 3 = at least 3
// The double ratio eps = (n_s2/n_no2)(deficit runs) / (n_s2/n_no2)(reference runs) needs no luminosity;
// wp2_runs.py groups the runs in time blocks.
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp2", kTRUE)' 'plots/all_report/wp2_runs.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>
#include "TSystem.h"
#include "defs.h"

namespace wp2r {
inline float Phi360(float phi) { return phi < 0.f ? phi + 360.f : phi; }
inline int Sec(float phi) {
    if (!(phi == phi)) return -1;
    phi = Phi360(phi);
    int k = (int)std::floor(phi / 60.);
    return ((k + 5) % 6 + 6) % 6;
}
struct C { double n[3][7] = {}, v[3][7] = {}; int day = 0; };   // [sample][tot, s2, no2, A12, A34, B012, B3]
}  // namespace wp2r

void wp2_runs(int sel = 0 /* 0 S0, 1 pre, 2 none */) {
    using namespace wp2r;
    scan::Data D;
    if (!scan::Load(D)) return;
    const int selId[3] = {alldefs::kSelS0, alldefs::kSelPre, alldefs::kSelNone};
    std::map<int, C> M;
    for (int smp = 0; smp < 3; ++smp) {
        const std::vector<scan::Row>& v = D.s[smp];
        int cur = -1, run = -1;
        double g[7] = {0, 0, 0, 0, 0, 0, 0};
        auto flush = [&]() {
            if (run >= 0) for (int k = 0; k < 7; ++k) M[run].v[smp][k] += g[k] * g[k];
            for (int k = 0; k < 7; ++k) g[k] = 0;
        };
        for (size_t i = 0; i < v.size(); ++i) {
            const scan::Row& r = v[i];
            if (r.grp != cur) { flush(); cur = r.grp; run = r.run; }
            if (!alldefs::PassSel(selId[sel], r, r.cor_ok == 1, r.dz_pi, r.mm4)) continue;
            const int s[4] = {Sec(r.pip_phi), Sec(r.pim_phi), Sec(r.ep_phi), Sec(r.em_phi)};
            const bool has2 = s[0] == 2 || s[1] == 2 || s[2] == 2 || s[3] == 2;
            C& c = M[r.run];
            c.day = r.file_id;
            c.n[smp][0] += 1; g[0] += 1;
            if (has2) { c.n[smp][1] += 1; g[1] += 1; } else { c.n[smp][2] += 1; g[2] += 1; }
            int use = 0;
            for (int k = 0; k < 4; ++k) use |= 1 << s[k];
            int kOther = 0;
            for (int q = 0; q < 6; ++q) if (q != 2 && (use & (1 << q))) ++kOther;
            const int cat = has2 ? (kOther >= 3 ? 6 : 5) : (kOther >= 3 ? 4 : 3);
            c.n[smp][cat] += 1; g[cat] += 1;
        }
        flush();
    }
    const char* tag[3] = {"S0", "pre", "none"};
    FILE* o = fopen(Form("/hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_runs_%s.csv", tag[sel]), "w");
    fprintf(o, "run,day");
    const char* sm[3] = {"us", "pp", "mm"};
    const char* kn[7] = {"tot", "s2", "no2", "A12", "A34", "B012", "B3"};
    for (int smp = 0; smp < 3; ++smp) for (int k = 0; k < 7; ++k) fprintf(o, ",%s_%s,%s_v%s", sm[smp], kn[k], sm[smp], kn[k]);
    fprintf(o, "\n");
    for (auto& x : M) {
        fprintf(o, "%d,%03d", x.first, x.second.day);
        for (int smp = 0; smp < 3; ++smp)
            for (int k = 0; k < 7; ++k) fprintf(o, ",%.0f,%.0f", x.second.n[smp][k], x.second.v[smp][k]);
        fprintf(o, "\n");
    }
    fclose(o);
    printf("wp2_runs: %zu runs written (%s)\n", M.size(), tag[sel]);
}
