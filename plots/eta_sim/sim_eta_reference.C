// ============================================================================
// sim_eta_reference.C  --  reference table of the simulated components with the eta set
// ============================================================================
// e+e- hypothesis (trees std_us, lep_us, eta_us), owned rows (std and lep in
// partition mode A, eta owned rows), weight w (PDG BR). For the selections
// "pre", "pre+OA>2", "slides" and "old code" (SimEngine presets) and three M4
// ranges (all, bump window 0.455-0.5225, peak window 0.5336-0.5542), per
// component:
//   rows (W > 0), sum w, N_eff (event groups), sigma [nb] = sum w / L_sim,
//   expected data rows = sum w * L_data / L_sim (L_data = 6080.5 nb^-1, UNCONFIRMED)
// Components: bump pi0 Dalitz, bump conversion, peak conversion, peak direct
// (eta -> pi+pi-e+e-, eta set), peak total, X17 for BR = 1e-5 (illustration,
// not in any sum), the signal sum, the physical background split by set (std,
// lep, eta; eta also its four-particle class 7 alone), the full simulated
// spectrum T1 + T2 (all sets) and its std + lep part (the first study), the
// ghost variant (a) and, as information, fake lepton and different sources.
// Rows with W = 0 (eta-set rows whose only forced decay is X17) are not
// counted in any component except the X17 illustration.
//
// Output: plots/eta_sim/tables/eta_engine_reference.txt
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_eta", kTRUE)' 'plots/eta_sim/sim_eta_reference.C+'
// ============================================================================
#include <cstdio>
#include <string>
#include <vector>

#include "TStopwatch.h"

#include "SimEngine.h"
#include "sim_ghostvar.h"

namespace seref {
using namespace simeng;
// drop rows with W = 0 (X17 products in the eta set) from a component (the X17 illustration keeps them)
inline Comp NZ(const Comp& c) {
    Comp o{c.name, [c](const Row& r) { return (r.w > 0 || r.w_db > 0) && c.f(r); }};
    o.wf = c.wf;
    return o;
}
}  // namespace seref

void sim_eta_reference() {
    using namespace simeng;
    using namespace seref;
    TStopwatch swAll;
    SimData d;
    if (!Load(d, kCompact, kMaskUS)) return;

    std::vector<Comp> comps;
    std::vector<std::string> label;
    auto add = [&](const Comp& c, const char* l) { comps.push_back(c); label.push_back(l); };
    add(NZ(CompBumpDalitz()), "bump, pi0 Dalitz (C2)");
    add(NZ(CompBumpConv()), "bump, conversion (C3)");
    add(NZ(CompPeakConv()), "peak, conversion (C4)");
    add(NZ(CompPeakDirect()), "peak, direct (C5, M(ee) not QED)");
    add(NZ(CompPeakTotal()), "peak total (C4 + C5)");
    add(CompX17(1e-5), "X17, BR 1e-5 (C6; illustration)");
    add(NZ({"signal", [](const Row& r) { return r.tmpl == 1 && !X17Row(r); }}), "signal sum T1 (C2-C5)");
    add(NZ(CompPhysBkgSets(kSetStd)), "phys. background, std (T2)");
    add(NZ(CompPhysBkgSets(kSetLep)), "phys. background, lep (T2)");
    add(NZ(CompPhysBkgSets(kSetEta)), "phys. background, eta (T2)");
    add(NZ({"eta C7", [](const Row& r) { return r.set == 3 && r.tmpl == 2 && r.cls == 7; }}), "  of which eta four-class 7");
    add(NZ(CompPhysBkg()), "phys. background, all sets");
    add(NZ(CompSim()), "T1+T2 full simulated spectrum");
    add(NZ(CompInSets(CompSim(), kSetStd | kSetLep)), "T1+T2, std + lep only (first study)");
    add(NZ(CompSimVar(1)), "variant (a): T1+T2+ghost true pairs");
    add(NZ(CompFake()), "info: fake lepton (T3), all sets");
    add(NZ(CompDiff()), "info: different sources (T4)");
    const int nc = (int)comps.size();

    const int presets[4] = {kPre, kPreOA2, kSlides, kOldCode};
    const int nr = 3;
    const double rlo[nr] = {-1e9, scan::kBumpLo, scan::kPeakLo};
    const double rhi[nr] = {1e9, scan::kBumpHi, scan::kPeakHi};
    const char* rname[nr] = {"all M4", "bump window 0.455-0.5225", "peak window 0.5336-0.5542"};

    const std::string outName = std::string(kTabDir) + "/eta_engine_reference.txt";
    FILE* o = fopen(outName.c_str(), "w");
    if (!o) { printf("cannot write %s\n", outName.c_str()); return; }
    auto both = [&](const std::string& s) { fputs(s.c_str(), o); fputs(s.c_str(), stdout); };
    both("# eta_engine_reference: simulated pi+pi-e+e- components (e+e- hypothesis) with the SMASH eta set\n");
    both("# sets: std (owned, mode A) + lep (owned, mode A) + eta (owned: >= 1 track of a forced eta decay); weight w (PDG BR)\n");
    both(Form("# sigma = sum w / L_sim, L_sim = %.2f nb^-1 for each set (the same 10^10 SMASH events); within the HADES acceptance\n", kLsim));
    both(Form("#   after the selection, no efficiency correction; data rows = sum w * L_data / L_sim, L_data = %.1f nb^-1 (UNCONFIRMED)\n", kLdata));
    both("# errors: event-group variance V = sum_g (sum w of the group)^2 (sets independent); N_eff = (sum w)^2 / V\n");
    both("# rows: rows with W > 0 (eta-set rows whose only forced decay is X17 have W = 0 and are not counted, except in the X17 line)\n");
    both("# CAUTION (review 2026-11-12, finding 2): the eta-set eta -> pi+pi-e+e- has a NON-QED generator M(ee) (accepted density largest at\n");
    both("#   threshold, close to 1/M^2; tables/eta_peak_mee_true.txt). Every 'peak direct' number below (yield, direct/conversion ratio, window\n");
    both("#   contents) is labelled 'generator M(ee) model not QED' and is not for physics until the PLUTO model is checked or reweighted.\n");
    both("# peak direct: W = 2 x BR(eta -> pi+pi-e+e-) = 2 x 2.68e-4 per row (q = 1/2); X17: W = 2 x 1e-5 per class-6 row,\n");
    both("#   INFORMATION ONLY (illustrative BR of the order of the BESIII limit 1.5e-5), never in a sum\n");
    both("# variant (a): T1 + T2 + ghost-flagged true one-source pairs (reconcile.md); fake lepton and different sources: information\n");
    both("# vertex: vz > -160 mm in every preset; slides = pre + oa_lab < 50, M(pipi) < 0.46, oa_eta > 140, OA > 2;\n");
    both("#   old code = pre + oa_lab < 50, M(pipi) < 0.42, oa_eta > 140\n");
    TStopwatch swFill;
    for (int ip = 0; ip < 4; ++ip) {
        const int p = presets[ip];
        for (int ir = 0; ir < nr; ++ir) {
            const double a = rlo[ir], b = rhi[ir];
            const Sel base = SelPreset(p);
            const Sel sel = ir == 0 ? base : SelAnd(base, [a, b](const Row& r) { return r.m4 >= a && r.m4 < b; });
            const std::vector<WSum> res = SumMany(d, sel, comps);   // default: std + lep + eta, mode A, w
            both(Form("\n## selection \"%s\", %s\n", PresetName(p), rname[ir]));
            both(Form("%-40s %10s %12s %11s %22s %24s\n", "component", "rows", "sum w", "N_eff", "sigma [nb]", "data rows (6080.5)"));
            for (int c = 0; c < nc; ++c) {
                const WSum& s = res[c];
                both(Form("%-40s %10.0f %12.6g %11.1f %10.4g +- %-9.3g %11.5g +- %-9.3g\n", label[c].c_str(), s.N, s.S, s.Neff(),
                          s.Xsec(), s.XsecErr(), s.DataRows(), std::sqrt(s.V) * kLdata / kLsim));
            }
            // shares in this range
            const WSum& sim = res[12];
            if (sim.S > 0)
                both(Form("# shares of T1+T2: peak direct %.4f, peak conversion %.4f, bump %.4f, phys. bkg eta set %.4f; direct / conversion %.3f\n",
                          res[3].S / sim.S, res[2].S / sim.S, (res[0].S + res[1].S) / sim.S, res[9].S / sim.S,
                          res[2].S > 0 ? res[3].S / res[2].S : 0.));
        }
    }
    // rules v2 (2026-11-13): 727.319 (rules v1: 727.299)
    both("\n# checks: 'pre', all M4, eta owned rows (W > 0 or not) must equal tables/eta_validate.txt (727.319 sum w);\n");
    both("#   T1+T2 std + lep only must equal tables/engine_reference.txt (698.9 nb in 'pre')\n");
    printf("sim_eta_reference: filling %.1f s, total %.1f s, table %s\n", swFill.RealTime(), swAll.RealTime(), outName.c_str());
    fclose(o);
}
