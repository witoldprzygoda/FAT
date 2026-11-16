// ============================================================================
// wp2_simfs.C  --  WP2: share of the simulated signal with a track in sector 2 (f_s of C_sec)
// ============================================================================
// Plain words: the simulation has no sector-2 deficit. For every simulated row of the nominal
// selection S0 (mode A, w x V3 = physics code 108, M4' = 0.99615 M4) we ask whether one of its four
// tracks lies in sector 2 (sector = (floor(phi/60) + 5) mod 6, phi brought to [0, 360)), and in which
// sector pattern (A_k = no sector-2 track, k sectors; B_k = sector-2 track(s) and k other sectors).
// f_s = weighted share of the rows with a sector-2 track, per signal class and window:
//   bump C2 + C3 (and each) in B; peak C4 + C5 (and each) in P; T12 (simulated sum) in B, P, CR-U and all M4.
// Also the share for every sector (the simulation is six-fold symmetric, so all should agree).
// Errors: event-group variance of the ratio (sum over groups of (w_cat - f w_all)^2 / W_all^2).
// Streaming read (no Load into memory): std_us, lep_us (sim_compact.root), eta_us (sim_compact_eta.root).
// Output: docs/latex/tables_all/wp2_simfs.txt (+ machine-readable lines "fs|...").
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp2", kTRUE)' \
//           'plots/all_report/wp2_simfs.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "defs.h"
#include "sim_eta_physics.h"

namespace wp2f {
inline float Phi360(float phi) { return phi < 0.f ? phi + 360.f : phi; }
inline int Sec(float phi) {
    if (!(phi == phi)) return -1;
    phi = Phi360(phi);
    int k = (int)std::floor(phi / 60.);
    return ((k + 5) % 6 + 6) % 6;
}
// groups of components
enum { gC23 = 0, gC2, gC3, gC45, gC4, gC5, gGeta, gT12, gBkg, kNG };
const char* kGName[kNG] = {"C2+C3", "C2", "C3", "C4+C5", "C4", "C5", "G_eta", "T12", "T12-C2..C5"};
enum { wB = 0, wP, wCRU, wAll, kNWw };
const char* kWName[kNWw] = {"B", "P", "CR-U", "all M4"};
// categories: 0 all, 1 uses s2, 2..7 uses sector 0..5, 8..11 A1..A4, 12..15 B0..B3, 16 one s2 track, 17 >=2 s2 tracks
const int kNCat = 18;
const char* kCName[kNCat] = {"all", "uses s2", "uses s0", "uses s1", "uses s2 (copy)", "uses s3", "uses s4", "uses s5",
                             "A1", "A2", "A3", "A4", "B0", "B1", "B2", "B3", "one s2 track", ">=2 s2 tracks"};
const int kN = kNG * kNWw * kNCat;
inline int Idx(int g, int w, int c) { return (g * kNWw + w) * kNCat + c; }
struct Acc {
    std::vector<double> W, Wg, Sc2, Sa2, Sca;   // totals; per-group buffer; sums of group products
    std::vector<int> touched;
    Acc() : W(kN, 0.), Wg(kN, 0.), Sc2(kN, 0.), Sa2(kN, 0.), Sca(kN, 0.) {}
    void Add(int i, double w) { if (Wg[i] == 0.) touched.push_back(i); Wg[i] += w; W[i] += w; }
    void Flush() {
        // the 'all' cell of each (g, w) is category 0
        for (int i : touched) {
            const int a = i - (i % kNCat);
            Sc2[i] += Wg[i] * Wg[i];
            Sca[i] += Wg[i] * Wg[a];
        }
        for (int i : touched) if (i % kNCat == 0) Sa2[i] += Wg[i] * Wg[i];
        for (int i : touched) Wg[i] = 0.;
        touched.clear();
    }
};
}  // namespace wp2f

void wp2_simfs() {
    using namespace wp2f;
    using namespace alldefs;
    TStopwatch sw;
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    Acc A;
    long long nRows = 0, nPass = 0, nBadPhi = 0;
    const char* files[3] = {simeng::kCompact, simeng::kCompact, simeng::kCompactEta};
    const char* trees[3] = {"std_us", "lep_us", "eta_us"};
    for (int t = 0; t < 3; ++t) {
        TFile* f = TFile::Open(files[t], "READ");
        TTree* tr = f ? (TTree*)f->Get(trees[t]) : nullptr;
        if (!tr) { printf("wp2_simfs: tree %s missing\n", trees[t]); return; }
        simeng::Row r;
        simeng::BindRow(tr, r, false);
        const Long64_t n = tr->GetEntries();
        int cur = -1;
        for (Long64_t i = 0; i < n; ++i) {
            tr->GetEntry(i);
            if (r.grp != cur) { A.Flush(); cur = r.grp; }
            ++nRows;
            if (r.ownA != 1) continue;
            const float dzpi = DzPi(r.pip_z, r.pim_z, r.vz);
            if (!PassSel(kSelS0, r, r.cor_ok == 1, dzpi, r.mm4)) continue;
            int cn[6];
            const int ncn = CompsOf(r, false, cn);
            if (ncn == 0) continue;
            const double w = r.w * simphys::RowFactor(r, kPhysNom);
            if (w == 0.) continue;
            ++nPass;
            const double m4p = M4sim(r.m4, r.m4_true);
            int wins[kNWw], nw = 0;
            if (m4p >= kBlo && m4p < kBhi) wins[nw++] = wB;
            if (m4p >= kPlo && m4p < kPhi) wins[nw++] = wP;
            if (m4p >= 0.565 && m4p < 0.700) wins[nw++] = wCRU;
            wins[nw++] = wAll;
            const float phi[4] = {r.pip_phi, r.pim_phi, r.ep_phi, r.em_phi};
            int use = 0, n2 = 0;
            bool bad = false;
            for (int k = 0; k < 4; ++k) {
                const int s = Sec(phi[k]);
                if (s < 0) { bad = true; break; }
                use |= 1 << s;
                if (s == 2) ++n2;
            }
            if (bad) { ++nBadPhi; continue; }
            int kOther = 0;
            for (int q = 0; q < 6; ++q) if (q != 2 && (use & (1 << q))) ++kOther;
            int cats[kNCat], nc = 0;
            cats[nc++] = 0;
            if (n2 > 0) cats[nc++] = 1;
            for (int s = 0; s < 6; ++s) if (use & (1 << s)) cats[nc++] = 2 + s;
            if (n2 == 0) cats[nc++] = 8 + (kOther - 1);
            else { cats[nc++] = 12 + kOther; cats[nc++] = n2 == 1 ? 16 : 17; }
            // groups of this row
            int gs[kNG], ng = 0;
            for (int q = 0; q < ncn; ++q) {
                switch (cn[q]) {
                    case kC2: gs[ng++] = gC23; gs[ng++] = gC2; gs[ng++] = gBkg; break;   // gBkg corrected below
                    case kC3: gs[ng++] = gC23; gs[ng++] = gC3; gs[ng++] = gBkg; break;
                    case kC4: gs[ng++] = gC45; gs[ng++] = gC4; gs[ng++] = gBkg; break;
                    case kC5: gs[ng++] = gC45; gs[ng++] = gC5; gs[ng++] = gBkg; break;
                    case kGeta: gs[ng++] = gGeta; break;
                    case kT12: gs[ng++] = gT12; break;
                    default: break;
                }
            }
            // gBkg = T12 minus the signal classes: a row that is a signal class gets no gBkg entry; a T2 row does
            bool sig = false;
            for (int q = 0; q < ncn; ++q) if (cn[q] >= kC2 && cn[q] <= kC5) sig = true;
            int gs2[kNG], ng2 = 0;
            for (int q = 0; q < ng; ++q) if (gs[q] != gBkg) gs2[ng2++] = gs[q];
            bool hasT12 = false;
            for (int q = 0; q < ncn; ++q) if (cn[q] == kT12) hasT12 = true;
            if (hasT12 && !sig) gs2[ng2++] = gBkg;
            for (int q = 0; q < ng2; ++q)
                for (int a = 0; a < nw; ++a)
                    for (int c = 0; c < nc; ++c) A.Add(Idx(gs2[q], wins[a], cats[c]), w);
        }
        A.Flush();
        tr->ResetBranchAddresses();
        f->Close();
        delete f;
        printf("wp2_simfs: %s done, %.1f s\n", trees[t], sw.RealTime());
        sw.Continue();
    }
    FILE* o = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp2_simfs.txt", "w");
    fprintf(o, "# wp2_simfs.C: share f of the simulated rows with a track in sector 2 (and other patterns), S0, SMASH std+lep+eta owned\n");
    fprintf(o, "# rows (mode A), w x V3 (108), M4' = 0.99615 M4; windows on M4'; errors = event-group variance of the ratio\n");
    fprintf(o, "# rows read %lld, passing S0 with a component %lld, NaN phi %lld\n", nRows, nPass, nBadPhi);
    fprintf(o, "# group   window   category        W [sum w]      share    +- err\n");
    for (int g = 0; g < kNG; ++g) {
        for (int w = 0; w < kNWw; ++w) {
            const int ia = Idx(g, w, 0);
            const double Wa = A.W[ia];
            if (Wa <= 0) continue;
            for (int c = 0; c < kNCat; ++c) {
                if (c == 4) continue;
                const int i = Idx(g, w, c);
                const double f = A.W[i] / Wa;
                const double var = (c == 0) ? 0. : (A.Sc2[i] - 2 * f * A.Sca[i] + f * f * A.Sa2[ia]) / (Wa * Wa);
                fprintf(o, "  %-10s %-7s %-16s %12.4f   %.4f +- %.4f\n", kGName[g], kWName[w], kCName[c], A.W[i], f, std::sqrt(std::max(0., var)));
                fprintf(o, "fs|%s|%s|%s,%.6f,%.6f,%.6f\n", kGName[g], kWName[w], kCName[c], f, std::sqrt(std::max(0., var)), A.W[i]);
            }
        }
    }
    fclose(o);
    printf("wp2_simfs: done in %.1f s; table docs/latex/tables_all/wp2_simfs.txt\n", sw.RealTime());
}
