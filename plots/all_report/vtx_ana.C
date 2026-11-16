// ==============================================================================================================
// vtx_ana.C -- copy (2026-11-14, integration of the event-vertex efficiency into report_all / report_exp) of
//   /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_verify/vv_ana.C
//   (study of the eta efficiency of the mandatory event-vertex cut vz > -500 mm; verifier ve_verify, own code).
// What it does: S (P_main) in V and Nv, model grid, free fits, like-sign null tests, ECAL-tagged T, bump edge step E.
// Changes against the original: file / function / include names (vv_* -> vtx_*), and the text outputs go to
//   /hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/vtx_*.txt instead of the session scratchpad. Numbers, cuts and fits are unchanged.
// The original outputs are copied to docs/latex/tables_all/vtx_*.txt; the notes to docs/latex/notes_all/ve_*_note.md.
// Run (from /home/przygoda/HADES/FAT/FAT; vtx_fill.C first, then the analysis macros, at most 6 processes at once):
//   B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic", kTRUE)'
//   root -l -b -q -e "$B" 'plots/all_report/vtx_ana.C+'
// ==============================================================================================================
// vtx_ana.C -- ve_verify: eta content of rows with (V) and without (Nv) an event vertex, own code
// root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic",kTRUE)' 'vtx_ana.C+'
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "vtx_lib.h"
using namespace vvl;

static FILE* gO = nullptr;
#define P(...) do { printf(__VA_ARGS__); if (gO) fprintf(gO, __VA_ARGS__); } while (0)

std::vector<R> Load(const char* pat, bool sim = false) {
  TChain ch("t"); ch.Add(pat);
  Int_t smp, day, grp, best, start, vcand, cor, ng; Float_t vz, m4, mee, mpp, oa, meeg, m5c, zx, rx, dca, zxe, rxe, z[4], r[4];
  ch.SetBranchStatus("*", 0);
  auto BI = [&](const char* n, Int_t* a) { ch.SetBranchStatus(n, 1); ch.SetBranchAddress(n, a); };
  auto BF = [&](const char* n, Float_t* a) { ch.SetBranchStatus(n, 1); ch.SetBranchAddress(n, a); };
  BI("smp", &smp); BI("day", &day); BI("grp", &grp); BI("best", &best); BI("start", &start); BI("vcand", &vcand);
  BI("cor", &cor); BI("ng", &ng);
  BF("vz", &vz); BF("m4", &m4); BF("mee", &mee); BF("mpp", &mpp); BF("oa", &oa); BF("meeg", &meeg); BF("m5c", &m5c);
  BF("zx", &zx); BF("rx", &rx); BF("dca", &dca); BF("zxe", &zxe); BF("rxe", &rxe); BF("z", z); BF("r", r);
  std::vector<R> v; v.reserve(ch.GetEntries());
  for (Long64_t i = 0; i < ch.GetEntries(); ++i) {
    ch.GetEntry(i);
    R q; q.file = ch.GetTreeNumber(); q.grp = grp; q.smp = smp; q.day = day; q.best = best; q.start = start; q.vcand = vcand;
    q.ng = ng; q.vz = vz; q.m4 = m4; q.mee = mee; q.mpp = mpp; q.oa = oa; q.meeg = meeg; q.m5c = m5c; q.zx = zx; q.rx = rx;
    q.dca = dca; q.zxe = zxe; q.rxe = rxe; q.cor = cor;
    float zs[4] = {z[0], z[1], z[2], z[3]}; std::sort(zs, zs + 4); q.zp = 0.5f * (zs[1] + zs[2]);
    int nw = 0; for (int k = 0; k < 4; ++k) if (z[k] > -170 && z[k] < -30 && std::fabs(r[k]) < 20) ++nw; q.nwin = nw;
    v.push_back(q);
  }
  return v;
}

struct C { std::string name; Sel s; };
int Period(int d) { return d <= 46 ? 0 : d <= 48 ? 1 : d <= 53 ? 2 : 3; }

void vtx_ana() {
  const char* outTxt = "/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/vtx_ana.txt";
  gO = fopen(outTxt, "w");
  std::vector<R> d = Load("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/rows/d_*.root");
  long long nt[3] = {0, 0, 0}, nv[3] = {0, 0, 0}, ntc[3] = {0, 0, 0}, nvc[3] = {0, 0, 0};
  for (const R& r : d) { nt[r.smp]++; if (r.vz > -500) nv[r.smp]++; if (r.cor) { ntc[r.smp]++; if (r.vz > -500) nvc[r.smp]++; } }
  P("# vtx_ana: rows read from the 29 PAT day files (isBest >= 0, trigbit == 8192)\n");
  for (int s = 0; s < 3; ++s)
    P("  smp %d: rows %lld, vz>-500 %lld (%.4f); cor_ok %lld, cor_ok & vz>-500 %lld (%.4f)\n", s, nt[s], nv[s], (double)nv[s] / nt[s],
      ntc[s], nvc[s], (double)nvc[s] / ntc[s]);

  // ---------------- classes ----------------

  std::vector<C> cl;
  auto V = [](const R& r) { return r.vz > -500; };
  auto Nv = [](const R& r) { return r.vz <= -500; };
  cl.push_back({"ALL", [](const R&) { return true; }});
  cl.push_back({"V", V});
  cl.push_back({"Nv", Nv});
  cl.push_back({"Vpre(vz>-160)", [](const R& r) { return r.vz > -160; }});
  cl.push_back({"V,-500<vz<=-160", [](const R& r) { return r.vz > -500 && r.vz <= -160; }});
  cl.push_back({"V,best1", [](const R& r) { return r.vz > -500 && r.best == 1; }});
  cl.push_back({"Nv,best1", [](const R& r) { return r.vz <= -500 && r.best == 1; }});
  const char* pn[4] = {"I", "IIA", "IIB", "III"};
  for (int p = 0; p < 4; ++p) {
    cl.push_back({Form("V,per%s", pn[p]), [p](const R& r) { return r.vz > -500 && Period(r.day) == p; }});
    cl.push_back({Form("Nv,per%s", pn[p]), [p](const R& r) { return r.vz <= -500 && Period(r.day) == p; }});
  }
  for (int st = 0; st < 2; ++st) {
    cl.push_back({Form("V,start%d", st), [st](const R& r) { return r.vz > -500 && r.start == st; }});
    cl.push_back({Form("Nv,start%d", st), [st](const R& r) { return r.vz <= -500 && r.start == st; }});
  }
  // pion crossing: on axis r<5, 5-8, >=8 mm (dca<10), no good crossing (dca>=10 or parallel)
  const char* rn[4] = {"rx<5", "5<=rx<8", "rx>=8", "dca>=10"};
  auto RxC = [](const R& r) { if (!(r.dca >= 0 && r.dca < 10)) return 3; return r.rx < 5 ? 0 : r.rx < 8 ? 1 : 2; };
  for (int k = 0; k < 4; ++k) {
    cl.push_back({Form("V,%s", rn[k]), [k, RxC](const R& r) { return r.vz > -500 && RxC(r) == k; }});
    cl.push_back({Form("Nv,%s", rn[k]), [k, RxC](const R& r) { return r.vz <= -500 && RxC(r) == k; }});
  }
  // median track z
  const char* zn[4] = {"zp<-300", "-300<=zp<-170", "-170<=zp<-30", "zp>=-30"};
  auto ZpC = [](const R& r) { return r.zp < -300 ? 0 : r.zp < -170 ? 1 : r.zp < -30 ? 2 : 3; };
  for (int k = 0; k < 4; ++k) {
    cl.push_back({Form("V,%s", zn[k]), [k, ZpC](const R& r) { return r.vz > -500 && ZpC(r) == k; }});
    cl.push_back({Form("Nv,%s", zn[k]), [k, ZpC](const R& r) { return r.vz <= -500 && ZpC(r) == k; }});
  }
  // own tracks in the finder window
  for (int k = 0; k < 3; ++k) {
    cl.push_back({Form("V,nwin%s", k == 2 ? ">=2" : k == 1 ? "=1" : "=0"), [k](const R& r) { return r.vz > -500 && std::min((int)r.nwin, 2) == k; }});
    cl.push_back({Form("Nv,nwin%s", k == 2 ? ">=2" : k == 1 ? "=1" : "=0"), [k](const R& r) { return r.vz <= -500 && std::min((int)r.nwin, 2) == k; }});
  }
  // M(ee)
  const char* mn[3] = {"mee<0.03", "0.03<=mee<0.10", "mee>=0.10"};
  auto MeC = [](const R& r) { return r.mee < 0.03 ? 0 : r.mee < 0.10 ? 1 : 2; };
  for (int k = 0; k < 3; ++k) {
    cl.push_back({Form("V,%s", mn[k]), [k, MeC](const R& r) { return r.vz > -500 && MeC(r) == k; }});
    cl.push_back({Form("Nv,%s", mn[k]), [k, MeC](const R& r) { return r.vz <= -500 && MeC(r) == k; }});
  }
  std::vector<Sel> sels; for (auto& c : cl) sels.push_back(c.s);
  std::vector<Spec> sp = Fill(d, sels);

  P("\n# A. P_main (frozen 0.5439 / 5.4 MeV, pol2, 0.5225-0.61, N - CB, CB 10 MeV interpolated), per class\n");
  P("  %-22s %9s %9s %7s %9s %8s %8s %6s %10s %7s %8s %8s %7s\n", "class", "N_win", "CB_win", "CB/N", "S", "stat", "stat_g", "chi2",
    "ndf", "p", "Spois*", "", "");
  std::vector<Fit> fits;
  for (size_t i = 0; i < cl.size(); ++i) {
    PCfg c; Fit f = Peak(sp[i], c); fits.push_back(f);
    double N, CB; vvl::WinSums(sp[i], 0.53364, 0.55416, N, CB);
    const double g = sp[i].winN > 0 ? std::sqrt(sp[i].winVarG / sp[i].winN) : 1;
    P("  %-22s %9.0f %9.0f %7.3f %9.0f %8.0f %8.0f %6.1f %4d %7.3f   (rows us/pp/mm %lld/%lld/%lld)\n", cl[i].name.c_str(), N, CB,
      N > 0 ? CB / N : 0, f.S(), f.eS(), f.eS() * g, f.chi2, f.ndf, f.prob, sp[i].nsel[0], sp[i].nsel[1], sp[i].nsel[2]);
  }
  auto idxOf = [&](const std::string& n) { for (size_t i = 0; i < cl.size(); ++i) if (cl[i].name == n) return (int)i; return -1; };
  auto gfac = [&](int i) { return sp[i].winN > 0 ? std::sqrt(sp[i].winVarG / sp[i].winN) : 1.0; };
  auto epsPair = [&](const std::string& a, const std::string& b) {
    const int ia = idxOf(a), ib = idxOf(b); double e;
    const double v = EpsV(fits[ia].S(), fits[ia].eS() * gfac(ia), fits[ib].S(), fits[ib].eS() * gfac(ib), e);
    P("  eps_S(%s | %s) = %.3f +- %.3f\n", a.c_str(), b.c_str(), v, e);
  };
  P("\n# B. eps_vtx from S = S(V)/(S(V)+S(Nv)), stat_g errors\n");
  epsPair("V", "Nv"); epsPair("V,best1", "Nv,best1");
  for (int p = 0; p < 4; ++p) epsPair(Form("V,per%s", pn[p]), Form("Nv,per%s", pn[p]));
  for (int k = 0; k < 4; ++k) epsPair(Form("V,%s", rn[k]), Form("Nv,%s", rn[k]));
  for (int k = 0; k < 3; ++k) epsPair(Form("V,%s", mn[k]), Form("Nv,%s", mn[k]));
  P("  sum check: S(V)+S(Nv) = %.0f vs S(ALL) = %.0f\n", fits[idxOf("V")].S() + fits[idxOf("Nv")].S(), fits[idxOf("ALL")].S());

  // ---------------- robustness of S(Nv) and S(V) ----------------
  const int iV = idxOf("V"), iN = idxOf("Nv");
  P("\n# C. model variations (S_V, S_Nv, eps)\n");
  P("  %-40s %8s %6s %8s %6s %7s %6s  chi2V/chi2N\n", "variation", "S_V", "e", "S_Nv", "e", "eps", "e");
  auto var = [&](const char* name, PCfg c) {
    Fit a = Peak(sp[iV], c), b = Peak(sp[iN], c); double e;
    const double v = EpsV(a.S(), a.eS() * gfac(iV), b.S(), b.eS() * gfac(iN), e);
    P("  %-40s %8.0f %6.0f %8.0f %6.0f %7.3f %6.3f  %.1f/%.1f (ndf %d)\n", name, a.S(), a.eS(), b.S(), b.eS(), v, e, a.chi2, b.chi2, a.ndf);
    return std::make_pair(a.S(), b.S());
  };
  std::vector<double> sNv, sV;
  for (int o = 1; o <= 4; ++o)
    for (double hi : {0.58, 0.59, 0.61, 0.63, 0.65})
      for (double lo : {0.5225, 0.5250, 0.5275})
        for (double ds : {-0.0005, 0.0, 0.0005}) {
          PCfg c; c.ord = o; c.hi = hi; c.lo = lo; c.sg = 0.0054 + ds;
          Fit a = Peak(sp[iV], c), b = Peak(sp[iN], c);
          if (a.prob > 0.05 && b.prob > 0.05) { sV.push_back(a.S()); sNv.push_back(b.S()); }
        }
  auto rms = [](const std::vector<double>& x, double& mean) { mean = 0; for (double v : x) mean += v; mean /= x.size();
    double s = 0; for (double v : x) s += (v - mean) * (v - mean); return std::sqrt(s / x.size()); };
  double mV, mN; const double rV = rms(sV, mV), rN = rms(sNv, mN);
  std::vector<double> ep; for (size_t i = 0; i < sV.size(); ++i) ep.push_back(sV[i] / (sV[i] + sNv[i]));
  double mE; const double rE = rms(ep, mE);
  P("  grid (order 1-4, hi 0.58-0.65, lo 0.5225-0.5275, width +-0.5 MeV; both p>0.05): %zu fits; S_V mean %.0f rms %.0f; S_Nv mean %.0f rms %.0f; eps mean %.3f rms %.3f (min %.3f max %.3f)\n",
    sV.size(), mV, rV, mN, rN, mE, rE, *std::min_element(ep.begin(), ep.end()), *std::max_element(ep.begin(), ep.end()));
  { PCfg c; var("nominal", c); }
  { PCfg c; c.sg = 0.0065; var("width 6.5 MeV", c); }
  { PCfg c; c.sg = 0.0075; var("width 7.5 MeV", c); }
  for (double k : {0.85, 0.92, 1.08, 1.15}) { PCfg c; c.k = k; var(Form("CB x %.2f", k), c); }
  { PCfg c; c.what = 4; Fit a = Peak(sp[iV], c), b = Peak(sp[iN], c);
    P("  free CB scale (fit N = S G + pol2 + k CB): V: S %.0f +- %.0f, k %.3f +- %.3f; Nv: S %.0f +- %.0f, k %.3f +- %.3f\n",
      a.S(), a.eS(), a.p[4], std::sqrt(a.cov(4, 4)), b.S(), b.eS(), b.p[4], std::sqrt(b.cov(4, 4))); }

  // free position and width
  P("\n# D. free position and width (profile chi2 grid, mean 0.534-0.554, sigma 3-12 MeV)\n");
  for (const char* n : {"ALL", "V", "Nv", "Vpre(vz>-160)", "Nv,perIII", "Nv,zp<-300", "Nv,-300<=zp<-170"}) {
    const int i = idxOf(n); Free f = FreeFit(sp[i], PCfg());
    P("  %-20s mean %.5f +- %.5f  sigma %.2f +- %.2f MeV  S %.0f +- %.0f  chi2 %.1f (frozen chi2 %.1f)%s\n", n, f.mu, f.emu, f.sg * 1e3,
      f.esg * 1e3, f.S, f.eS, f.chi2, fits[i].chi2, f.edgeHit ? "  [at grid edge]" : "");
  }

  // like-sign null tests
  P("\n# E. like-sign null tests: frozen Gaussian + pol2 on N++, N--, N++ + N-- (Poisson variance)\n");
  for (const char* n : {"ALL", "V", "Nv"}) {
    const int i = idxOf(n);
    PCfg c1; c1.what = 1; PCfg c2; c2.what = 2; PCfg c3; c3.what = 3;
    Fit a = Peak(sp[i], c1), b = Peak(sp[i], c2), s = Peak(sp[i], c3);
    P("  %-6s N++: %.0f +- %.0f   N--: %.0f +- %.0f   sum: %.0f +- %.0f   (2 sqrt-weighted CB-equivalent ~ %.0f)\n", n, a.S(), a.eS(), b.S(),
      b.eS(), s.S(), s.eS(), 0.0);
  }

  // look-elsewhere: frozen-width fits at displaced positions, same relative range
  P("\n# F. displaced-position scan (frozen width 5.4 MeV, range m-21.4..m+66.1 MeV, pol2, N - CB): pulls S/stat\n");
  for (const char* n : {"V", "Nv"}) {
    const int i = idxOf(n); std::vector<double> pulls; std::string line;
    for (double m = 0.575; m <= 0.8651; m += 0.0075) {
      PCfg c; c.mu = m; c.lo = m - 0.0214; c.hi = m + 0.0661; Fit f = Peak(sp[i], c);
      if (!f.ok) continue; pulls.push_back(f.S() / f.eS()); line += Form(" %.1f", f.S() / f.eS());
    }
    double mp; const double rp = rms(pulls, mp);
    int n2 = 0; for (double p : pulls) if (std::fabs(p) > 2) ++n2;
    P("  %-3s positions 0.575-0.865 step 7.5 MeV: %zu fits, mean pull %.2f, rms %.2f, |pull|>2: %d\n    pulls:%s\n", n, pulls.size(), mp, rp, n2, line.c_str());
    // the eta position itself
    PCfg c; Fit f = Peak(sp[i], c); P("    at eta: pull %.2f\n", f.S() / f.eS());
  }

  // ---------------- ECAL-tagged eta T in m5c ----------------
  P("\n# G. ECAL-tagged eta T: tag 0.119 <= M(ee g) < 0.173, M(ee) < 0.10; fit N - CB in m5c 0.49-0.62, Gaussian 0.5463 / 6.2 MeV + pol2\n");
  {
    auto tag = [](const R& r) { return r.ng > 0 && r.meeg >= 0.119 && r.meeg < 0.173 && r.mee < 0.10; };
    std::vector<std::string> tn = {"ALL", "V", "Nv", "Vpre(vz>-160)", "Nv,perIII", "V,perIII", "Nv,start1", "Nv,zp>=-170"};
    std::vector<Sel> ts;
    ts.push_back([tag](const R& r) { return tag(r); });
    ts.push_back([tag](const R& r) { return tag(r) && r.vz > -500; });
    ts.push_back([tag](const R& r) { return tag(r) && r.vz <= -500; });
    ts.push_back([tag](const R& r) { return tag(r) && r.vz > -160; });
    ts.push_back([tag](const R& r) { return tag(r) && r.vz <= -500 && Period(r.day) == 3; });
    ts.push_back([tag](const R& r) { return tag(r) && r.vz > -500 && Period(r.day) == 3; });
    ts.push_back([tag](const R& r) { return tag(r) && r.vz <= -500 && r.start == 1; });
    ts.push_back([tag](const R& r) { return tag(r) && r.vz <= -500 && r.zp >= -170; });
    std::vector<Spec> st = Fill(d, ts, [](const R& r) { return r.m5c; }, 0.5463 - 1.9 * 0.0062, 0.5463 + 1.9 * 0.0062);
    std::vector<Fit> tf;
    for (size_t i = 0; i < ts.size(); ++i) {
      PCfg c; c.mu = 0.5463; c.sg = 0.0062; c.lo = 0.49; c.hi = 0.62; Fit f = Peak(st[i], c); tf.push_back(f);
      Free ff = FreeFit(st[i], c);
      P("  %-16s T = %7.0f +- %4.0f  chi2 %.1f/%d  (rows us %lld)  free: mean %.4f sigma %.1f MeV S %.0f\n", tn[i].c_str(), f.S(), f.eS(), f.chi2,
        f.ndf, st[i].nsel[0], ff.mu, ff.sg * 1e3, ff.S);
    }
    double e; const double v = EpsV(tf[1].S(), tf[1].eS(), tf[2].S(), tf[2].eS(), e);
    P("  eps_T = %.3f +- %.3f\n", v, e);
    // pi0 tag rate per (N - CB) of M(ee) < 0.10 rows in V and Nv
    std::vector<Sel> ms = {[](const R& r) { return r.vz > -500 && r.mee < 0.10; }, [](const R& r) { return r.vz <= -500 && r.mee < 0.10; },
                           [tag](const R& r) { return r.vz > -500 && tag(r); }, [tag](const R& r) { return r.vz <= -500 && tag(r); }};
    std::vector<Spec> sm = Fill(d, ms);
    double a[4];
    for (int k = 0; k < 4; ++k) { double N = 0, CB = 0; for (int i = 0; i < kNb; ++i) { N += sm[k].n[0][i]; CB += sm[k].cb[i]; } a[k] = N - CB;
      P("  M4 0.35-0.95, %s: N - CB = %.0f (N %.0f)\n", k == 0 ? "V, mee<0.1" : k == 1 ? "Nv, mee<0.1" : k == 2 ? "V, tagged" : "Nv, tagged", N - CB, N); }
    P("  tag rate per (N - CB): V %.4f, Nv %.4f\n", a[2] / a[0], a[3] / a[1]);
  }

  // ---------------- bump edge step E ----------------
  P("\n# H. bump edge step E (4 M(pipi) slices, edge fixed, pol2 + [eta] + h step, half range 45/45/45/35 MeV)\n");
  {
    const double sl[5][2] = {{0.30, 0.34}, {0.34, 0.37}, {0.37, 0.39}, {0.39, 0.40}, {0.42, 0.46}};
    std::vector<Sel> es;
    for (int k = 0; k < 5; ++k)
      for (int c = 0; c < 3; ++c) {
        const double lo = sl[k][0], hi = sl[k][1];
        es.push_back([lo, hi, c](const R& r) { return r.mpp >= lo && r.mpp < hi && (c == 0 ? r.vz > -500 : c == 1 ? r.vz <= -500 : r.vz > -160); });
      }
    std::vector<Spec> se = Fill(d, es);
    double E[3] = {0, 0, 0}, E2[3] = {0, 0, 0};
    const char* cn[3] = {"V", "Nv", "Vpre"};
    // model variations: order 2/3, half range 35/45/55, smearing 4.4/5.4/6.4
    std::vector<std::vector<double>> varE(3);
    for (int k = 0; k < 5; ++k) {
      for (int c = 0; c < 3; ++c) {
        ECfg ec; ec.mpLo = sl[k][0]; ec.mpHi = sl[k][1]; if (k == 3) ec.half = 0.035; if (k == 4) { ec.useEdge = true; ec.edge = 0.5150; }
        bool we; Fit f = Edge(se[k * 3 + c], ec, we);
        P("  slice %.2f-%.2f %-5s h = %6.0f +- %4.0f  chi2 %.1f/%d %s\n", sl[k][0], sl[k][1], cn[c], f.S(), f.eS(), f.chi2, f.ndf, we ? "(eta in model)" : "");
        if (k < 4) { E[c] += f.S(); E2[c] += f.eS() * f.eS(); }
      }
    }
    for (int c = 0; c < 3; ++c) P("  E(%s) = %.0f +- %.0f\n", cn[c], E[c], std::sqrt(E2[c]));
    double e; const double v = EpsV(E[0], std::sqrt(E2[0]), E[1], std::sqrt(E2[1]), e);
    P("  eps_E = %.3f +- %.3f\n", v, e);
    // variations applied to both classes
    std::vector<double> epsv;
    for (int ord = 2; ord <= 3; ++ord) for (int ih = 0; ih < 3; ++ih) for (int is = -1; is <= 1; ++is) {
      double Ev[2] = {0, 0};
      for (int k = 0; k < 4; ++k) for (int c = 0; c < 2; ++c) {
        ECfg ec; ec.mpLo = sl[k][0]; ec.mpHi = sl[k][1]; ec.ord = ord; ec.half = 0.035 + 0.010 * ih; ec.s = 0.0054 + 0.001 * is;
        bool we; Fit f = Edge(se[k * 3 + c], ec, we); Ev[c] += f.S(); }
      epsv.push_back(Ev[0] / (Ev[0] + Ev[1]));
      P("    var ord %d half %.3f s %.4f: E_V %.0f E_Nv %.0f eps %.3f\n", ord, 0.035 + 0.010 * ih, 0.0054 + 0.001 * is, Ev[0], Ev[1], Ev[0] / (Ev[0] + Ev[1]));
    }
    double me; const double re = rms(epsv, me); P("  eps_E over 18 variations: mean %.3f rms %.3f\n", me, re);
  }
  fclose(gO);
}
