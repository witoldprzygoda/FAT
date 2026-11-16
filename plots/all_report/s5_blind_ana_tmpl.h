// ============================================================================
// s5_blind_ana_tmpl.h  --  half-A templates of the blind test, as used by the
// blind analysis (label blind_ana). Reads ONLY blind/templates/templates.root.
// ============================================================================
// Atomic (mutually exclusive) components per selection, M4 in 1 MeV bins at
// L_sim: direct5, conv4, bump2, bump3, Geta, Gpi0D, Gconv, Gvec, Goth.
// Sums: peak45 = direct5 + conv4, bump23 = bump2 + bump3,
//       bkgT2 = Geta + Gpi0D + Gconv + Gvec + Goth, sumT12 = all.
// Every number is scaled by `scale` (0.5 for the eta toys at L_sim / 2).
// Replicas: each atomic bin fluctuated by a Gaussian with its event-group error
// (template MC statistics), sums rebuilt from the fluctuated atoms.
// ============================================================================
#ifndef S5_BLIND_ANA_TMPL_H
#define S5_BLIND_ANA_TMPL_H

#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"

namespace bla {

const char* const kTemplFile = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/templates/templates.root";
const char* const kAtoms[9] = {"direct5", "conv4", "bump2", "bump3", "Geta", "Gpi0D", "Gconv", "Gvec", "Goth"};
enum { kDirect5 = 0, kConv4, kBump2, kBump3, kGeta, kGpi0D, kGconv, kGvec, kGoth, kNAtom };
const char* const kSel[4] = {"PEAK", "DIRECT", "BUMP", "X17SEL"};
const int kNSel = 4;

// study windows (simulation-centred, fixed for all agents)
const double kPeakLo = 0.5333, kPeakHi = 0.5587;
const double kBumpLo = 0.4576, kBumpHi = 0.5255;
const double kPeakMean = 0.5460, kPeakSigma = 0.0067;   // simulated peak (frozen)

// M4 template arrays on the 1 MeV grid 0.30-1.20 (900 bins), at L_sim
struct M4Templ {
    std::vector<double> c[kNAtom], e2[kNAtom];   // content and variance per 1 MeV bin
    double regPeak[kNAtom], regBump[kNAtom];     // exact window sums (reg_ histograms)
    double regPeakE[kNAtom], regBumpE[kNAtom];
};

inline bool LoadM4(TFile& f, int sel, M4Templ& t) {
    for (int a = 0; a < kNAtom; ++a) {
        TH1D* h = (TH1D*)f.Get(Form("m4_%s_%s", kSel[sel], kAtoms[a]));
        TH1D* r = (TH1D*)f.Get(Form("reg_%s_%s", kSel[sel], kAtoms[a]));
        if (!h || !r) { printf("missing template %s %s\n", kSel[sel], kAtoms[a]); return false; }
        t.c[a].assign(900, 0.); t.e2[a].assign(900, 0.);
        for (int b = 1; b <= 900; ++b) { t.c[a][b - 1] = h->GetBinContent(b); t.e2[a][b - 1] = h->GetBinError(b) * h->GetBinError(b); }
        t.regBump[a] = r->GetBinContent(2); t.regBumpE[a] = r->GetBinError(2);
        t.regPeak[a] = r->GetBinContent(4); t.regPeakE[a] = r->GetBinError(4);
    }
    return true;
}

// sum of atoms (mask of atom indices) on a grid of width `rb` x 1 MeV starting at 0.30, scaled
inline std::vector<double> Sum(const M4Templ& t, const std::vector<int>& atoms, int rb, double scale, std::vector<double>* var = nullptr,
                               const std::vector<std::vector<double>>* repl = nullptr) {
    const int nb = 900 / rb;
    std::vector<double> s(nb, 0.);
    if (var) var->assign(nb, 0.);
    for (int a : atoms)
        for (int b = 0; b < 900; ++b) {
            const double v = repl ? (*repl)[a][b] : t.c[a][b];
            s[b / rb] += scale * v;
            if (var) (*var)[b / rb] += scale * scale * t.e2[a][b];
        }
    return s;
}
// content between x1 and x2 of a sum of atoms (1 MeV bins, fractional edge bins), scaled
inline double SumWin(const M4Templ& t, const std::vector<int>& atoms, double x1, double x2, double scale, double* err = nullptr,
                     const std::vector<std::vector<double>>* repl = nullptr) {
    double s = 0, e2 = 0;
    for (int a : atoms)
        for (int b = 0; b < 900; ++b) {
            const double lo = 0.30 + 0.001 * b, hi = lo + 0.001;
            const double l = std::max(lo, x1), h = std::min(hi, x2);
            if (h <= l) continue;
            const double fr = (h - l) / 0.001;
            s += scale * fr * (repl ? (*repl)[a][b] : t.c[a][b]);
            e2 += scale * scale * fr * fr * t.e2[a][b];
        }
    if (err) *err = std::sqrt(e2);
    return s;
}
// one replica of the atoms: Gaussian fluctuation of each bin by its error
inline std::vector<std::vector<double>> Replica(const M4Templ& t, TRandom3& rnd) {
    std::vector<std::vector<double>> r(kNAtom);
    for (int a = 0; a < kNAtom; ++a) {
        r[a].resize(900);
        for (int b = 0; b < 900; ++b) r[a][b] = t.c[a][b] + (t.e2[a][b] > 0 ? rnd.Gaus(0, std::sqrt(t.e2[a][b])) : 0.);
    }
    return r;
}

inline std::vector<int> AtPeak() { return {kDirect5, kConv4}; }
inline std::vector<int> AtBump() { return {kBump2, kBump3}; }
inline std::vector<int> AtBkg() { return {kGeta, kGpi0D, kGconv, kGvec, kGoth}; }
inline std::vector<int> AtBkgNoEta() { return {kGpi0D, kGconv, kGvec, kGoth}; }
inline std::vector<int> AtAll() { return {kDirect5, kConv4, kBump2, kBump3, kGeta, kGpi0D, kGconv, kGvec, kGoth}; }

} // namespace bla
#endif
