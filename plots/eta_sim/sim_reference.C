// ============================================================================
// sim_reference.C  --  reference table of the simulated pi+pi-e+e- rows
// ============================================================================
// For the selections "pre" and "pre+OA>2" (SimEngine presets), per template
// class and per four-particle class of the e+e- hypothesis (trees std_us and
// lep_us, owned rows of both sets added):
//   rows, sum w, N_eff (event groups), cross section in nb (sum w / 221.05)
//   and the expected rows in data for L_data = 6080.5 nb^-1 (unconfirmed),
// in partition mode A and mode B, with the weights w (PDG) and w_db; the
// signal classes also without the rows with a ghost lepton; the split of
// the default variant into the std and lep sets.
//
// Output: plots/eta_sim/tables/engine_reference.txt
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_engine", kTRUE)' 'plots/eta_sim/sim_reference.C+'
// ============================================================================
#include <cstdio>
#include <string>
#include <vector>

#include "TStopwatch.h"

#include "SimEngine.h"

void sim_reference() {
    using namespace simeng;
    TStopwatch swAll;
    SimData d;
    if (!Load(d, kCompact, (1 << kStdUS) | (1 << kLepUS))) return;
    TStopwatch swFill;

    // components: templates, sum, four-particle classes, signal classes without ghost leptons
    std::vector<Comp> comps;
    for (int t = 1; t <= 4; ++t) comps.push_back(CompTmpl(t));
    comps.push_back(CompSim());
    comps.push_back(CompAny());
    for (int c = 0; c <= 8; ++c) comps.push_back(CompCls(c));
    for (int c = 2; c <= 4; ++c) comps.push_back(CompSigCls(c, false));
    for (int c = 2; c <= 4; ++c) comps.push_back(CompSigCls(c, true));
    std::vector<std::string> label;
    for (int t = 1; t <= 4; ++t) label.push_back(Form("T%d %s", t, TmplName(t)));
    label.push_back("T1+T2 simulated spectrum");
    label.push_back("all owned rows (T1-T4)");
    for (int c = 0; c <= 8; ++c) label.push_back(Form("C%d %s", c, ClsName(c)));
    for (int c = 2; c <= 4; ++c) label.push_back(Form("signal C%d %s", c, ClsName(c)));
    for (int c = 2; c <= 4; ++c) label.push_back(Form("signal C%d, no ghost lepton", c));

    const int presets[2] = {kPre, kPreOA2};
    // variants: (mode, db)
    const int nv = 4;
    const int vMode[nv] = {0, 0, 1, 1};
    const bool vDb[nv] = {false, true, false, true};
    const char* vName[nv] = {"A, w", "A, w_db", "B, w", "B, w_db"};

    const std::string outName = std::string(kTabDir) + "/engine_reference.txt";
    FILE* o = fopen(outName.c_str(), "w");
    if (!o) { printf("cannot write %s\n", outName.c_str()); return; }
    auto both = [&](const char* s) { fputs(s, o); fputs(s, stdout); };

    both("# engine_reference: simulated pi+pi-e+e- rows (e+e- hypothesis), std + lep owned rows added\n");
    both(Form("# sigma = sum w / L_sim, L_sim = %.2f nb^-1; expected data rows = sum w * L_data / L_sim, L_data = %.1f nb^-1 (UNCONFIRMED)\n",
              kLsim, kLdata));
    both("# errors: event-group variance V = sum_g (sum w of the group)^2; N_eff = (sum w)^2 / V; no efficiency correction\n");
    both("# mode A: lep provides rows with a forced pi0/eta Dalitz (or vector meson) product; mode B: std keeps them, lep only rho0/omega/phi\n");
    both("# w: PDG branching ratios (default); w_db: converter-database BR, no natural-decay factors\n");

    for (int ip = 0; ip < 2; ++ip) {
        const int p = presets[ip];
        const Sel sel = SelPreset(p);
        std::vector<std::vector<WSum>> res(nv);
        for (int v = 0; v < nv; ++v) {
            SimOpt opt; opt.mode = vMode[v]; opt.db = vDb[v];
            res[v] = SumMany(d, sel, comps, opt);
        }
        // per set, default variant
        SimOpt oStd; oStd.sets = 1;
        SimOpt oLep; oLep.sets = 2;
        const std::vector<WSum> rStd = SumMany(d, sel, comps, oStd), rLep = SumMany(d, sel, comps, oLep);

        both(Form("\n## selection \"%s\"\n", PresetName(p)));
        both("### default: mode A, w (PDG)\n");
        both(Form("%-40s %10s %12s %12s %16s %20s\n", "class", "rows", "sum w", "N_eff", "sigma [nb]", "data rows (6080.5)"));
        for (size_t c = 0; c < comps.size(); ++c) {
            const WSum& s = res[0][c];
            both(Form("%-40s %10.0f %12.6g %12.1f %8.4g +- %-6.2g %11.4g +- %-7.2g\n", label[c].c_str(), s.N, s.S, s.Neff(),
                      s.Xsec(), s.XsecErr(), s.DataRows(), std::sqrt(s.V) * kLdata / kLsim));
        }
        both("### variants: sigma [nb] +- stat (rows)\n");
        both(Form("%-40s", "class"));
        for (int v = 0; v < nv; ++v) both(Form(" %28s", vName[v]));
        both("\n");
        for (size_t c = 0; c < comps.size(); ++c) {
            both(Form("%-40s", label[c].c_str()));
            for (int v = 0; v < nv; ++v)
                both(Form("   %9.4g +- %-7.2g (%7.0f)", res[v][c].Xsec(), res[v][c].XsecErr(), res[v][c].N));
            both("\n");
        }
        both("### default variant split by set: sigma [nb] +- stat (rows, N_eff)\n");
        both(Form("%-40s %36s %36s\n", "class", "std owned", "lep owned"));
        for (size_t c = 0; c < comps.size(); ++c)
            both(Form("%-40s   %9.4g +- %-7.2g (%7.0f, %9.1f)   %9.4g +- %-7.2g (%7.0f, %9.1f)\n", label[c].c_str(),
                      rStd[c].Xsec(), rStd[c].XsecErr(), rStd[c].N, rStd[c].Neff(), rLep[c].Xsec(), rLep[c].XsecErr(),
                      rLep[c].N, rLep[c].Neff()));
    }
    // rules v2 (2026-11-13): 112014 + 72804 = 184818 (rules v1: 112030 + 72804 = 184834)
    both(Form("\n# check (\"pre\", all owned rows): mode A sum w must equal 112014 (std) + 72804 (lep) = 184818 of the validated production\n"));
    printf("sim_reference: filling %.1f s, total %.1f s, table %s\n", swFill.RealTime(), swAll.RealTime(), outName.c_str());
    fclose(o);
}
