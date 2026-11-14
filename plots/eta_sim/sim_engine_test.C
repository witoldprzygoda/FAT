// ============================================================================
// sim_engine_test.C  --  self-test of SimEngine.h against independent numbers
// ============================================================================
//  1. sum w and N_eff (event groups) of all owned rows after "pre", per set,
//     against the validated partition tables (part_std_us.txt, part_lep_us.txt)
//  2. precision: the same sum accumulated in float and in a TH1F (the trap)
//  3. binned spectrum: sum over bins + under/overflow = total; variance >= sum w^2
//  4. efficiency: EfficiencyVs total = ratio of SumMany sums; error against the
//     uncorrelated formula of the spec
//  5. data path: DataSpectrum without vz cut (cor_ok only) in the ScanEngine
//     binning reproduces the data reference (peak window N 86011, CB 11488)
//  6. the presets on data: rows in the peak window per preset
// Output: plots/eta_sim/tables/engine_selftest.txt
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_engine", kTRUE)' 'plots/eta_sim/sim_engine_test.C+'
// ============================================================================
#include <cstdio>
#include <string>
#include <vector>

#include "TH1F.h"
#include "TStopwatch.h"

#include "SimEngine.h"

void sim_engine_test() {
    using namespace simeng;
    TStopwatch sw;
    FILE* o = fopen((std::string(kTabDir) + "/engine_selftest.txt").c_str(), "w");
    auto P = [&](const std::string& s) { fputs(s.c_str(), stdout); if (o) fputs(s.c_str(), o); };
    int nfail = 0;
    auto check = [&](bool ok, const std::string& what) { P(Form("  [%s] %s\n", ok ? "PASS" : "FAIL", what.c_str())); if (!ok) nfail++; };

    SimData d;
    if (!Load(d, kCompact, (1 << kStdUS) | (1 << kLepUS))) return;
    const Sel pre = SelPreset(kPre);

    // ---- 1 ----
    P("1. owned rows after \"pre\" against part_*_us.txt (mode A, w)\n");
    SimOpt oS; oS.sets = 1;
    SimOpt oL; oL.sets = 2;
    const WSum s1 = SumMany(d, pre, std::vector<Comp>(1, CompAny()), oS)[0];
    const WSum s2 = SumMany(d, pre, std::vector<Comp>(1, CompAny()), oL)[0];
    // references: rules v2 (2026-11-13), part_{std,lep}_us.txt of simrules.h v2 and tables/engine_validate_summary.txt
    // (rules v1: std sum w 112030, N_eff 93690.1; lep N_eff 5723165.7; mode B 193309 + 24.760 = 193334)
    P(Form("  std: rows %.0f, sum w %.8g, N_eff %.1f (reference 115414, 112014, 93784.0)\n", s1.N, s1.S, s1.Neff()));
    P(Form("  lep: rows %.0f, sum w %.8g, N_eff %.1f (reference 7082339, 72804, 5723161.1)\n", s2.N, s2.S, s2.Neff()));
    check(s1.N == 115414 && std::fabs(s1.S - 112014) < 1 && std::fabs(s1.Neff() - 93784.0) < 0.1, "std owned rows, sum w, N_eff(group)");
    check(s2.N == 7082339 && std::fabs(s2.S - 72804) < 1 && std::fabs(s2.Neff() - 5723161.1) < 0.2, "lep owned rows, sum w, N_eff(group)");
    SimOpt oB; oB.mode = 1;
    const WSum sB = SumMany(d, pre, std::vector<Comp>(1, CompAny()), oB)[0];
    P(Form("  mode B, std + lep: sum w %.8g (reference 193776 + 24.760 = 193800.8)\n", sB.S));
    check(std::fabs(sB.S - 193800.8) < 2, "mode B total");

    // ---- 2 ----
    P("2. precision of the accumulation (lep owned rows, pre)\n");
    float fsum = 0;
    TH1F hf("hf", "", 1, -1e9, 1e9);
    TH1D hd("hd", "", 1, -1e9, 1e9);
    for (const Row& r : d.t[kLepUS])
        if (r.ownA == 1 && pre(r)) { fsum += (float)r.w; hf.Fill(0., r.w); hd.Fill(0., r.w); }
    P(Form("  double (engine) %.8g ; float variable %.8g ; TH1F %.8g ; TH1D %.8g\n", s2.S, fsum, hf.GetBinContent(1), hd.GetBinContent(1)));
    check(std::fabs(hd.GetBinContent(1) - s2.S) < 1e-6 * s2.S, "TH1D equals the double sum");

    // ---- 3 ----
    P("3. binned spectrum consistency (signal + physical background, M4 0.3-1.1, 5 MeV)\n");
    Binning b; b.nbin = 160; b.lo = 0.3; b.hi = 1.1;
    const WSpec sp = Fill(d, pre, CompSim(), SimOpt(), b);
    double ss = 0, s2b = 0;
    for (int i = 0; i < sp.nbin; ++i) { ss += sp.s[i]; s2b += sp.v[i]; }
    P(Form("  sum of bins %.8g = S %.8g ; all rows %.8g (outside the range %.4g) ; sum V_bin %.6g, V(range) %.6g, V(all) %.6g\n",
           ss, sp.S, sp.Sall, sp.Sall - sp.S, s2b, sp.V, sp.Vall));
    check(std::fabs(ss - sp.S) < 1e-9 * sp.S && sp.Vall >= sp.V * (1 - 1e-12), "bin sums and variances");

    // ---- 4 ----
    P("4. efficiency of OA > 2 for the bump Dalitz signal and the physical background (pre)\n");
    const Sel pass = [](const Row& r) { return r.oa > 2.f; };
    for (const Comp& c : {CompBumpDalitz(), CompBumpConv(), CompPeakConv(), CompPhysBkg()}) {
        const EffRes e = Efficiency(d, pre, pass, c);
        const WSum a = SumMany(d, SelAnd(pre, pass), std::vector<Comp>(1, c))[0];
        const WSum f = SumMany(d, SelAnd(pre, [](const Row& r) { return !(r.oa > 2.f); }), std::vector<Comp>(1, c))[0];
        const EffRes u = EffFromSums(a.S, a.V, f.S, f.V);
        P(Form("  %-28s eps %.5f +- %.5f (group) ; from sums %.5f +- %.5f (no pass/fail correlation)\n", c.name.c_str(), e.eps, e.err, u.eps, u.err));
        check(std::fabs(e.eps - u.eps) < 1e-9, std::string("efficiency value, ") + c.name);
    }

    // ---- 5 ----
    P("5. data path against the data reference (no vz cut, cor_ok only; ScanEngine binning)\n");
    scan::Data dd;
    if (!scan::Load(dd)) return;
    const DSpec ds = DataSpectrum(dd, scan::SelNone(), scan::kNbin, scan::kMlo, scan::kMhi);
    const scan::WinSum& w = ds.sp.win[0];
    P(Form("  reb %d ; peak window N %.0f, CB %.0f (reference 86011, 11488)\n", ds.reb, w.n[0], w.cb));
    check(ds.reb == 4 && w.n[0] == 86011 && std::fabs(w.cb - 11488) < 1, "data N and CB in the peak window");
    TH1D* hs = ds.Hist("sig", "t5", 2);
    double integ = 0;
    for (int i = 1; i <= hs->GetNbinsX(); ++i) integ += hs->GetBinContent(i) * hs->GetBinWidth(i);
    double rows = 0;
    for (int i = 0; i < ds.sp.nbin; ++i) rows += ds.sp.sig[i];
    P(Form("  integral of dsigma/dM (N - CB) %.6g nb = rows %.6g / L_data %.1f = %.6g nb\n", integ, rows, kLdata, rows / kLdata));
    check(std::fabs(integ - rows / kLdata) < 1e-9 * std::fabs(integ), "data unit conversion");

    // ---- 6 ----
    P("6. presets on data: rows in the peak window 0.5336-0.5542 (N, CB)\n");
    std::vector<scan::Sel> sels;
    for (int p = 0; p < kNPreset; ++p) sels.push_back(DataSelPreset(p));
    scan::FillOpt fo;
    const std::vector<scan::Spectra> v = scan::FillMany(dd, sels, fo);
    for (int p = 0; p < kNPreset; ++p)
        P(Form("  %-26s N %7.0f  CB %7.0f  rows (any M4) %9.0f\n", PresetName(p), v[p].win[0].n[0], v[p].win[0].cb, v[p].nsel[0]));

    P(Form("\nself-test: %s, %.1f s\n", nfail ? "FAILURES" : "all checks passed", sw.RealTime()));
    if (o) fclose(o);
}
