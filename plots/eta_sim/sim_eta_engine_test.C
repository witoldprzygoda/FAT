// ============================================================================
// sim_eta_engine_test.C  --  self-test of the eta-set extension of SimEngine.h
// ============================================================================
//  1. backward compatibility: with only std_us + lep_us loaded, the default
//     SimOpt (sets = 7) gives the numbers of the first study (engine_selftest,
//     rules v2: std 115414 rows / 112014, lep 7082339 / 72804) and OptLabel says "std + lep"
//  2. eta trees: loaded owned rows, "pre" owned rows and sum w against
//     tables/eta_validate.txt (sim_eta_checks.C, read from scan_nt directly;
//     the "pre" and per-class / per-template references are read from that table)
//  3. per class and template of the eta set (rows, sum w) against eta_validate
//  4. default sum = std + lep + eta; variances add (separate trees)
//  5. X17: CompX17(BR) = 2 BR per class-6 row, linear in BR; never in CompSim
//  6. peak total = peak conversion + peak direct (sums and variances)
//  7. binned spectrum of the eta components: bin sums = SumMany; TH1D = double
//  8. efficiency with a weight override (X17) = ratio of SumMany sums
// Output: plots/eta_sim/tables/eta_engine_selftest.txt
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_eta", kTRUE)' 'plots/eta_sim/sim_eta_engine_test.C+'
// ============================================================================
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "TH1D.h"
#include "TStopwatch.h"

#include "SimEngine.h"

// reference numbers of the eta set: block [eta pippimepem] of tables/eta_validate.txt (sim_eta_checks.C, read from
// scan_nt directly): owned rows with vz > -160 and their sum w ('pre'), rows and sum w per four-particle class and
// per template. Missing values stay negative (the checks then fail).
struct EtaValRef {
    double preN = -1, preS = -1;
    std::map<int, std::pair<double, double> > cls, tmpl;
};
EtaValRef ReadEtaValidate() {
    EtaValRef r;
    std::ifstream in((std::string(simeng::kTabDir) + "/eta_validate.txt").c_str());
    std::string l;
    bool blk = false;
    while (std::getline(in, l)) {
        if (!l.empty() && l[0] == '[') { blk = l.compare(0, 16, "[eta pippimepem]") == 0; continue; }
        if (!blk) continue;
        const size_t p = l.find("owned with vz > -160 ");
        if (p != std::string::npos) sscanf(l.c_str() + p + 21, "%lf, sum w %lf", &r.preN, &r.preS);
        int k = 0;
        double n = 0, s = 0;
        if (sscanf(l.c_str(), " class %d: %lf %lf", &k, &n, &s) == 3) r.cls[k] = std::make_pair(n, s);
        if (sscanf(l.c_str(), " T%d: %lf %lf", &k, &n, &s) == 3) r.tmpl[k] = std::make_pair(n, s);
    }
    if (r.preS < 0 || r.cls.size() != 6 || r.tmpl.size() != 4)
        printf("sim_eta_engine_test: tables/eta_validate.txt incomplete ('pre' %s, %zu classes, %zu templates)\n", r.preS < 0 ? "missing" : "found",
               r.cls.size(), r.tmpl.size());
    return r;
}

void sim_eta_engine_test() {
    using namespace simeng;
    TStopwatch sw;
    FILE* o = fopen((std::string(kTabDir) + "/eta_engine_selftest.txt").c_str(), "w");
    auto P = [&](const std::string& s) { fputs(s.c_str(), stdout); if (o) fputs(s.c_str(), o); };
    int nfail = 0, npass = 0;
    auto check = [&](bool ok, const std::string& what) {
        P(Form("  [%s] %s\n", ok ? "PASS" : "FAIL", what.c_str()));
        if (ok) npass++; else nfail++;
    };
    const Sel pre = SelPreset(kPre);
    const std::vector<Comp> any(1, CompAny());

    // ---- 1 ----
    P("1. backward compatibility: only std_us + lep_us loaded, default SimOpt (sets = 7)\n");
    {
        SimData d0;
        if (!Load(d0, kCompact, (1 << kStdUS) | (1 << kLepUS))) return;
        SimOpt def;
        SimOpt o3; o3.sets = kSetStd | kSetLep;
        const WSum a = SumMany(d0, pre, any, def)[0], b = SumMany(d0, pre, any, o3)[0];
        P(Form("  default: rows %.0f, sum w %.8g, V %.8g ; sets = 3: rows %.0f, sum w %.8g, V %.8g ; label '%s'\n", a.N, a.S, a.V, b.N,
               b.S, b.V, OptLabel(def).c_str()));
        check(a.N == b.N && a.S == b.S && a.V == b.V, "default sets = 7 equals sets = 3 when no eta tree is loaded");
        // rules v2 (2026-11-13): 112014 + 72804 (engine_selftest, engine_validate_summary); rules v1: 112030 + 72804
        check(a.N == 115414 + 7082339 && std::fabs(a.S - (112014 + 72804)) < 2, "std + lep owned rows and sum w of the first study");
        check(OptLabel(def) == "std + lep, mode A, w (PDG)", "OptLabel without eta trees");
        check(d0.sets() == (kSetStd | kSetLep), "SimData::sets");
    }

    SimData d;
    if (!Load(d, kCompact, kMaskUS)) return;

    // ---- 2 ----
    P("2. eta trees against tables/eta_validate.txt\n");
    SimOpt oE; oE.sets = kSetEta;
    const EtaValRef ref = ReadEtaValidate();
    {
        P(Form("  eta_us loaded (owned A or B, all vz) %zu (reference 3188536)\n", d.t[kEtaUS].size()));
        check(d.t[kEtaUS].size() == 3188536, "eta_us owned rows loaded");
        const WSum e = SumMany(d, pre, any, oE)[0];
        P(Form("  'pre' owned: rows %.0f, sum w %.8g, N_eff %.1f (reference %.0f, %.8g)\n", e.N, e.S, e.Neff(), ref.preN, ref.preS));
        check(e.N == ref.preN && std::fabs(e.S - ref.preS) < 2e-4, "eta 'pre' owned rows and sum w");
        SimOpt oB = oE; oB.mode = 1;
        const WSum eb = SumMany(d, pre, any, oB)[0];
        check(eb.N == e.N && eb.S == e.S, "eta set owns the same rows in modes A and B");
        bool allSet3 = true;
        for (const Row& r : d.t[kEtaUS]) if (r.set != 3) { allSet3 = false; break; }
        check(allSet3, "set column == 3 in eta_us");
    }

    // ---- 3 ----
    P("3. eta set per four-particle class and template ('pre', owned)\n");
    {
        const int refCls[6] = {0, 1, 5, 6, 7, 8};
        double refN[6], refS[6];   // from tables/eta_validate.txt (rules v1 hard-coded values replaced on 2026-11-13)
        for (int k = 0; k < 6; ++k) {
            const auto it = ref.cls.find(refCls[k]);
            refN[k] = it == ref.cls.end() ? -1 : it->second.first;
            refS[k] = it == ref.cls.end() ? -1 : it->second.second;
        }
        std::vector<Comp> cc;
        for (int k = 0; k < 6; ++k) cc.push_back(CompCls(refCls[k]));
        for (int t = 1; t <= 4; ++t) cc.push_back(CompTmpl(t));
        const std::vector<WSum> r = SumMany(d, pre, cc, oE);
        bool ok = true;
        for (int k = 0; k < 6; ++k) {
            P(Form("  class %d: rows %9.0f (ref %9.0f), sum w %.8g (ref %.8g)\n", refCls[k], r[k].N, refN[k], r[k].S, refS[k]));
            if (r[k].N != refN[k] || std::fabs(r[k].S - refS[k]) > 1e-4 * std::max(1., refS[k])) ok = false;
        }
        check(ok, "rows and sum w per four-particle class");
        double refTN[4], refTS[4];
        for (int t = 0; t < 4; ++t) {
            const auto it = ref.tmpl.find(t + 1);
            refTN[t] = it == ref.tmpl.end() ? -1 : it->second.first;
            refTS[t] = it == ref.tmpl.end() ? -1 : it->second.second;
        }
        ok = true;
        for (int t = 0; t < 4; ++t) {
            P(Form("  T%d: rows %9.0f (ref %9.0f), sum w %.8g (ref %.8g)\n", t + 1, r[6 + t].N, refTN[t], r[6 + t].S, refTS[t]));
            if (r[6 + t].N != refTN[t] || std::fabs(r[6 + t].S - refTS[t]) > 1e-4 * refTS[t]) ok = false;
        }
        check(ok, "rows and sum w per template");
        const WSum dir = SumMany(d, pre, std::vector<Comp>(1, CompPeakDirect()), oE)[0];
        P(Form("  peak direct: sum w / rows = %.6g (2 x 2.68e-4 = 5.36e-4)\n", dir.S / dir.N));
        check(std::fabs(dir.S / dir.N - 5.36e-4) < 1e-9, "peak-direct weight per row = 2 BR");
    }

    // ---- 4 ----
    P("4. default sum = std + lep + eta\n");
    {
        SimOpt oS; oS.sets = kSetStd;
        SimOpt oL; oL.sets = kSetLep;
        const WSum all = SumMany(d, pre, any)[0], s = SumMany(d, pre, any, oS)[0], l = SumMany(d, pre, any, oL)[0],
                   e = SumMany(d, pre, any, oE)[0];
        P(Form("  all %.10g = std %.10g + lep %.10g + eta %.10g = %.10g ; V %.8g = %.8g ; label '%s'\n", all.S, s.S, l.S, e.S,
               s.S + l.S + e.S, all.V, s.V + l.V + e.V, OptLabel(SimOpt()).c_str()));
        check(std::fabs(all.S - (s.S + l.S + e.S)) < 1e-9 * all.S && std::fabs(all.V - (s.V + l.V + e.V)) < 1e-9 * all.V &&
                  all.N == s.N + l.N + e.N,
              "sum and variance of the three sets");
        check(OptLabel(SimOpt()) == "std + lep + eta, mode A, w (PDG)", "OptLabel with the eta trees");
    }

    // ---- 5 ----
    P("5. X17 illustration\n");
    {
        const std::vector<WSum> x = SumMany(d, pre, {CompX17(1e-5), CompX17(1.5e-5), CompCls(6)});
        P(Form("  CompX17(1e-5): rows %.0f, sum w %.8g (2e-5 x rows = %.8g); CompX17(1.5e-5) %.8g ; class-6 stored sum w %.3g\n",
               x[0].N, x[0].S, 2e-5 * x[0].N, x[1].S, x[2].S));
        // relative tolerance: 787 659 double additions round at the 1e-10 level
        check(x[0].N == 787659 && std::fabs(x[0].S / (2e-5 * x[0].N) - 1) < 1e-10 && std::fabs(x[1].S / x[0].S - 1.5) < 1e-10,
              "X17 weight 2 BR per class-6 row, linear in BR");
        check(x[2].S == 0, "stored X17 weights are 0 (BR(X17) = 0 in the production)");
        bool noX17 = true;
        const Comp cs = CompSim();
        for (const Row& r : d.t[kEtaUS]) if (r.cls == 6 && r.tmpl == 1 && cs.f(r)) { noX17 = false; break; }
        check(noX17, "CompSim never contains an X17 signal row");
    }

    // ---- 6 ----
    P("6. peak total = conversion + direct\n");
    {
        const std::vector<WSum> p = SumMany(d, pre, {CompPeakConv(), CompPeakDirect(), CompPeakTotal()});
        P(Form("  conv %.8g (V %.6g) + direct %.8g (V %.6g) = %.8g ; total %.8g (V %.6g)\n", p[0].S, p[0].V, p[1].S, p[1].V,
               p[0].S + p[1].S, p[2].S, p[2].V));
        check(std::fabs(p[2].S - p[0].S - p[1].S) < 1e-9 * p[2].S && std::fabs(p[2].V - p[0].V - p[1].V) < 1e-9 * p[2].V,
              "peak total sums and variances");
    }

    // ---- 7 ----
    P("7. binned spectra of the eta components (M4 0.3-1.1, 5 MeV)\n");
    {
        Binning b; b.nbin = 160; b.lo = 0.3; b.hi = 1.1;
        const std::vector<Comp> cc = {CompPeakDirect(), CompX17(1e-5), CompEtaBkg()};
        const std::vector<WSpec> sp = FillMany(d, pre, cc, SimOpt(), b);
        const std::vector<WSum> sm = SumMany(d, pre, cc);
        bool ok = true;
        for (size_t c = 0; c < cc.size(); ++c) {
            double ss = 0;
            for (int i = 0; i < sp[c].nbin; ++i) ss += sp[c].s[i];
            TH1D* h = sp[c].Hist(Form("t7_%zu", c), 0);
            P(Form("  %-40s bins %.8g, S %.8g, all rows %.8g (SumMany %.8g), TH1D integral %.8g\n", cc[c].name.c_str(), ss, sp[c].S,
                   sp[c].Sall, sm[c].S, h->Integral()));
            if (std::fabs(ss - sp[c].S) > 1e-9 * sp[c].S || std::fabs(sp[c].Sall - sm[c].S) > 1e-9 * sm[c].S ||
                std::fabs(h->Integral() - ss) > 1e-9 * ss)
                ok = false;
            delete h;
        }
        check(ok, "bin sums, totals and TH1D content");
    }

    // ---- 8 ----
    P("8. efficiency of OA > 2 with the weight override (X17) and for the direct peak\n");
    {
        const Sel pass = [](const Row& r) { return r.oa > 2.f; };
        bool ok = true;
        for (const Comp& c : {CompX17(1e-5), CompPeakDirect()}) {
            const EffRes e = Efficiency(d, pre, pass, c);
            const WSum a = SumMany(d, SelAnd(pre, pass), std::vector<Comp>(1, c))[0];
            const WSum all = SumMany(d, pre, std::vector<Comp>(1, c))[0];
            P(Form("  %-40s eps %.5f +- %.5f ; sums %.5f\n", c.name.c_str(), e.eps, e.err, a.S / all.S));
            if (std::fabs(e.eps - a.S / all.S) > 1e-9) ok = false;
        }
        check(ok, "efficiency equals the ratio of sums");
    }

    P(Form("\neta self-test: %d PASS, %d FAIL, %.1f s\n", npass, nfail, sw.RealTime()));
    if (o) fclose(o);
}
