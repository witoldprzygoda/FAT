// smp_report_cnt.C -- COPY of plots/exp_report/smp_report.C made by figure agent A (integ4, 2026-11-15) with
// the figure strings 'rows' -> 'counts' (analyst's terminology) and the outputs of the copy redirected:
// figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/smp_report
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/mk_cnt.py; the physics code is unchanged.
// ============================================================================
// smp_report.C  --  tables of the data sample, the luminosity and the full cut
//                   flow; figures of the vertex cut and of the cut flow (data)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Light macro: it only combines the outputs of
//   plots/exp_report/smp_pat_counts.C  (.../report_exp/smp/pat/pat_<day>.txt, .root, runs_<day>.txt)
//   plots/exp_report/smp_flow.C        (.../report_exp/smp/flow.txt)
// with the START-free luminosity tables (.../report_exp/smp/lumi_inputs/*.csv, revision N).
//
// What it does
//   1. data sample table: 29 day files (days, periods, file sizes, PAT rows of the
//      three hypothesis trees, EpEm_ID rows, runs of the luminosity table, runs
//      present in the data, L of all runs and of the analysed runs);
//   2. luminosity table per period: L6_SF of all runs with its errors, L of the
//      analysed runs (sum of L6_SF over the runs with EpEm_ID rows), the PT3 rules
//      (x 64/63; PT3/PT2 live fraction LF3/LF2 x 64/63 per day, L-weighted), L x
//      live fraction, the LGAD-START-condition luminosity of the published chain,
//      the four no_hf runs (L = 0) and their rows;
//   3. full cut flow: PAT rows -> isBest >= 0 -> vz > -500 mm -> PT3 (= FAT rows)
//      -> cor_ok -> vz > -160 mm (pre) -> pions in the vertex fit -> dz_pi < 40 mm,
//      and the named selections (each from pre), for e+e-, e+e+, e-e-: rows and
//      event groups for isBest >= 0, rows for isBest == 1; e+e- rows and CB/N in
//      the windows, without and with OA(e+e-) > 2 deg;
//   4. per-day rates against the PT3/PT2 live fraction: does dividing by
//      LF3/LF2 x 64/63 make the rates flatter? (correlation coefficients)
//   5. figures: smp_vertex_cut (eVertReco_z of the PAT rows and the fraction kept
//      by vz > -500 per day), smp_cutflow (rows and events kept per step).
// Selections: see smp_flow.C (pre = cor_ok && isBest >= 0 && vz > -160 mm, ...).
//
// Outputs: docs/latex/tables_exp/smp_samples.txt/.tex, smp_lumi.txt/.tex,
//          smp_cutflow.txt/.tex, smp_cutflow_windows.tex
//          docs/latex/figs_exp/smp_vertex_cut.pdf/.png, smp_cutflow.pdf/.png
//
// Run (from the repository root /home/przygoda/HADES/FAT/FAT; after smp_pat_counts.C and smp_flow.C; a few seconds):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_smp", kTRUE)' 'plots/exp_report/smp_report.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include "smp_common_cnt.h"

namespace smpr {

FILE* gOut = nullptr;
void P(const char* fmt, ...) {
    va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
    if (gOut) { va_start(a, fmt); vfprintf(gOut, fmt, a); va_end(a); }
}

// ---- PAT counts of one day ----
struct PatDay {
    int day = 0; long long size = 0, epem = 0;
    long long tree[3] = {0, 0, 0};
    long long rows[3][4], best[3][4], ev[3][4];
    long long vzNoVtx[3] = {0, 0, 0}, vzFail[3] = {0, 0, 0}, isBestNeg[3] = {0, 0, 0}, grpBestNot1[3] = {0, 0, 0};
    std::map<int, long long> trigBefore[3], trigAfter[3];
    int runsEpEm = 0, runsUS = 0, runsUS3 = 0;
    PatDay() { for (int k = 0; k < 3; ++k) for (int s = 0; s < 4; ++s) rows[k][s] = best[k][s] = ev[k][s] = 0; }
};
int Kof(const std::string& t) { return t == "us" ? 0 : (t == "pp" ? 1 : (t == "mm" ? 2 : -1)); }
bool ReadPat(int day, PatDay& p) {
    std::ifstream in(Form("%s/pat/pat_%03d.txt", smp::kWork, day));
    if (!in) return false;
    p.day = day;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream is(line);
        std::string key; is >> key;
        if (key == "FILE") { int d; is >> d >> p.size; }
        else if (key == "TREE") {
            std::string t; long long n; is >> t >> n;
            if (t == "EpEm_ID") p.epem = n;
            else if (t == "PipPimEpEm_ID") p.tree[0] = n;
            else if (t == "PipPimEpEp_ID") p.tree[1] = n;
            else if (t == "PipPimEmEm_ID") p.tree[2] = n;
        } else if (key == "STEP") {
            std::string t, nm; int s; long long a, b, c; is >> t >> s >> nm >> a >> b >> c;
            const int k = Kof(t); p.rows[k][s] = a; p.best[k][s] = b; p.ev[k][s] = c;
        } else if (key == "CHECK") {
            std::string t, w; long long a, b, c, d; is >> t >> w >> a >> w >> b >> w >> c >> w >> d;
            const int k = Kof(t); p.isBestNeg[k] = a; p.grpBestNot1[k] = b; p.vzNoVtx[k] = c; p.vzFail[k] = d;
        } else if (key == "TRIG") {
            std::string t, w; int tb; long long n; is >> t >> w >> tb >> n;
            const int k = Kof(t);
            if (w == "before") p.trigBefore[k][tb] += n; else p.trigAfter[k][tb] += n;
        } else if (key == "RUNS") {
            int d; std::string w; is >> d >> w >> p.runsEpEm >> w >> p.runsUS >> w >> p.runsUS3;
        }
    }
    return true;
}

// ---- flow.txt ----
struct CF { double rows = 0, ev = 0, w[4] = {0, 0, 0, 0}, cb[4] = {-1, -1, -1, -1}; std::string key; };
struct DayRate { double n0 = 0, v0 = 0, n1 = 0, R = 0, eR = 0, Rb = 0, eRb = 0; };

} // namespace smpr

void smp_report_cnt()
{
    using namespace smpr;
    smp::Style();
    gSystem->mkdir(smp::kWorkOut, kTRUE); gSystem->mkdir(smp::kTabOut, kTRUE); gSystem->mkdir(smp::kCntScratch, kTRUE);
    gOut = fopen(Form("%s/report.txt", smp::kWorkOut), "w");
    P("# smp_report.C (%s)\n", TDatime().AsString());
    const std::vector<int> days = smp::Days();
    const char* hn[3] = {"e+e-", "e+e+", "e-e-"};

    // ---------------- inputs ----------------
    std::map<int, PatDay> pat;
    for (int d : days) { PatDay p; if (!ReadPat(d, p)) { P("missing pat_%03d.txt\n", d); return; } pat[d] = p; }
    std::map<int, smp::DayLumi> dl = smp::ReadDayLumi();
    std::vector<smp::PerLumi> pl = smp::ReadPeriodLumi();
    std::map<int, smp::RunLumi> rl = smp::ReadRunLumi();
    std::map<int, smp::RunRows> rp = smp::ReadRunsPresent();
    std::map<int, double> La = smp::AnalysedLumiPerDay(false), Las = smp::AnalysedLumiPerDay(true);

    // flow.txt
    std::map<std::string, CF> cf;              // key "v i k"
    std::map<int, std::vector<DayRate>> dr;    // day -> [k]
    std::vector<std::string> stepKey(13);
    std::map<int, double> nohfRun;
    double nohf[2][3] = {{0, 0, 0}, {0, 0, 0}};
    {
        std::ifstream in(Form("%s/flow.txt", smp::kWork));
        if (!in) { P("flow.txt missing: run smp_flow.C first\n"); return; }
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream is(line);
            std::string key; is >> key;
            if (key == "CF") {
                int v, i, k; CF c; is >> v >> i >> c.key >> k >> c.rows >> c.ev >> c.w[0] >> c.w[1] >> c.w[2] >> c.w[3] >> c.cb[0] >> c.cb[1] >> c.cb[2] >> c.cb[3];
                cf[Form("%d %d %d", v, i, k)] = c;
                if (i < 13) stepKey[i] = c.key;
            } else if (key == "DAY") {
                int d, k; DayRate r; is >> d >> k >> r.n0 >> r.v0 >> r.n1 >> r.R >> r.eR >> r.Rb >> r.eRb;
                if (!dr.count(d)) dr[d].resize(3);
                dr[d][k] = r;
            } else if (key == "NOHF") {
                std::string w; is >> w;
                if (w == "rows") {
                    // NOHF rows FAT us a pp b mm c ; pre us d pp e mm f
                    std::string x; is >> x >> x >> nohf[0][0] >> x >> nohf[0][1] >> x >> nohf[0][2] >> x >> x >> x >> nohf[1][0] >> x >> nohf[1][1] >> x >> nohf[1][2];
                }
            } else if (key == "NOHFRUN") {
                int r; double n; is >> r >> n; nohfRun[r] = n;
            }
        }
    }
    auto C = [&](int v, int i, int k) -> const CF& { return cf[Form("%d %d %d", v, i, k)]; };

    // ---------------- 1. data sample table ----------------
    std::map<int, int> runsTab, runsNoHF, runsPresent, runsUS;
    for (const auto& p : rl) { runsTab[p.second.day]++; if (p.second.noHF) runsNoHF[p.second.day]++; }
    for (const auto& p : rp) {
        if (p.second.epem > 0) runsPresent[p.second.day]++;
        if (p.second.us3 > 0) runsUS[p.second.day]++;
    }
    {
        FILE* o = fopen(Form("%s/smp_samples.txt", smp::kTabOut), "w");
        fprintf(o, "# smp_samples.txt -- written by plots/exp_report/smp_report.C (%s)\n", TDatime().AsString());
        fprintf(o, "# FAT input: %s (29 day files lep_<day>.root, merged per day from the PAT outputs)\n", smp::kPatList);
        fprintf(o, "# PAT: pat3 tag prod-feb22-betanew-2026-10-28 (ana 0129bf3, scripts e54e6bd), Hydra 2.7.6, ROOT 6.36.02, Debian 12;\n");
        fprintf(o, "#      PID on beta_new with the polygons pp45_pid_cuts_exp_betanew_v1 (3 sigma p, pi+-, 5 sigma e+-), rejected hypotheses not written\n");
        fprintf(o, "# PAT rows = entries of the trees (one row per stored combination, all isBest >= 0); runs: luminosity table (lumi_coverage_perrun.csv),\n");
        fprintf(o, "#   present = runs with EpEm_ID rows in the day file, us = runs with FAT e+e- rows (after vz > -500, PT3)\n");
        fprintf(o, "# L6 = START-free luminosity of all runs of the day; L_an = sum of L6_SF over the runs present [nb^-1]\n");
        fprintf(o, "%-4s %-4s %8s | %10s %10s %10s %11s | %6s %5s %7s %6s | %9s %9s\n", "day", "per", "size[GB]", "PipPimEpEm", "PipPimEpEp", "PipPimEmEm",
                "EpEm_ID", "runs", "no_hf", "present", "us", "L6", "L_an");
        long long tsize = 0, tt[3] = {0, 0, 0}, te = 0; int tr = 0, tn = 0, tp = 0, tu = 0; double tL = 0, tLa = 0;
        for (int d : days) {
            const PatDay& p = pat[d];
            fprintf(o, "%03d  %-4s %8.2f | %10lld %10lld %10lld %11lld | %6d %5d %7d %6d | %9.3f %9.3f\n", d, smp::PeriodName(smp::PeriodOf(d)),
                    p.size / 1e9, p.tree[0], p.tree[1], p.tree[2], p.epem, runsTab[d], runsNoHF[d], runsPresent[d], runsUS[d], dl[d].L6, La[d]);
            tsize += p.size; for (int k = 0; k < 3; ++k) tt[k] += p.tree[k]; te += p.epem;
            tr += runsTab[d]; tn += runsNoHF[d]; tp += runsPresent[d]; tu += runsUS[d]; tL += dl[d].L6; tLa += La[d];
        }
        fprintf(o, "%-9s %8.2f | %10lld %10lld %10lld %11lld | %6d %5d %7d %6d | %9.3f %9.3f\n", "total", tsize / 1e9, tt[0], tt[1], tt[2], te, tr, tn, tp, tu, tL, tLa);
        // per period
        fprintf(o, "\n# per period\n");
        for (int q = 0; q < smp::kNPer; ++q) {
            long long s = 0, a[3] = {0, 0, 0}; int nr = 0, np = 0; double L = 0, Lx = 0; int nd = 0;
            for (int d : days) if (smp::PeriodOf(d) == q) {
                nd++; s += pat[d].size; for (int k = 0; k < 3; ++k) a[k] += pat[d].tree[k]; nr += runsTab[d]; np += runsPresent[d]; L += dl[d].L6; Lx += La[d];
            }
            fprintf(o, "%-4s %-8s %2d days %8.2f GB | %10lld %10lld %10lld | runs %5d present %5d | L6 %9.3f L_an %9.3f\n", smp::PeriodName(q),
                    smp::PeriodDays(q), nd, s / 1e9, a[0], a[1], a[2], nr, np, L, Lx);
        }
        fclose(o);
        P("samples: %lld bytes (%.1f GB), PAT rows us %lld pp %lld mm %lld, EpEm_ID %lld, runs %d (no_hf %d), present %d, with e+e- FAT rows %d, L6 %.3f, L_an %.3f\n",
          tsize, tsize / 1e9, tt[0], tt[1], tt[2], te, tr, tn, tp, tu, tL, tLa);

        // LaTeX
        FILE* x = fopen(Form("%s/smp_samples.tex", smp::kTabOut), "w");
        fprintf(x, "%% smp_samples.tex -- written by plots/exp_report/smp_report.C (numbers: tables_exp/smp_samples.txt)\n");
        fprintf(x, "\\begin{table}[H]\n\\centering\\footnotesize\n\\setlength{\\tabcolsep}{3.5pt}\n");
        fprintf(x, "\\caption{Experimental data sample: the 29 day files of the beam time (feb22, p+p, $T=4.53$~GeV) read by FAT "
                   "(PAT production tag \\texttt{prod-feb22-betanew-2026-10-28}, Hydra 2.7.6, PID on $\\beta_\\mathrm{new}$ with the polygons "
                   "\\texttt{betanew\\_v1}). PAT rows = stored particle combinations of the three hypothesis trees before any FAT cut. "
                   "Runs: in the luminosity table / present in the data (runs with $e^+e^-$ pair rows). "
                   "$L$: START-free luminosity (six sectors) of all runs and of the analysed runs, in nb$^{-1}$.}\n");
        fprintf(x, "\\label{tab:smp_samples}\n\\begin{tabular}{llr rrr rr rr}\n\\toprule\n");
        fprintf(x, " & & & \\multicolumn{3}{c}{PAT rows $\\pi^+\\pi^-\\ell\\ell$} & \\multicolumn{2}{c}{Runs} & \\multicolumn{2}{c}{$L$ [nb$^{-1}$]} \\\\\n");
        fprintf(x, "\\cmidrule(lr){4-6}\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}\n");
        fprintf(x, "Day & Period & GB & $e^+e^-$ & $e^+e^+$ & $e^-e^-$ & all & present & all & analysed \\\\\n\\midrule\n");
        for (int d : days) {
            const PatDay& p = pat[d];
            fprintf(x, "%03d & %s & %.1f & %lld & %lld & %lld & %d & %d & %.1f & %.1f \\\\\n", d, smp::PeriodName(smp::PeriodOf(d)), p.size / 1e9,
                    p.tree[0], p.tree[1], p.tree[2], runsTab[d], runsPresent[d], dl[d].L6, La[d]);
        }
        fprintf(x, "\\midrule\nTotal & & %.1f & %lld & %lld & %lld & %d & %d & %.1f & %.1f \\\\\n", tsize / 1e9, tt[0], tt[1], tt[2], tr, tp, tL, tLa);
        fprintf(x, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        fclose(x);
    }

    // ---------------- 2. luminosity ----------------
    {
        // per period: analysed L, live fraction, LGAD
        double Lan[4] = {0, 0, 0, 0}, Lans[4] = {0, 0, 0, 0}, Lall[4] = {0, 0, 0, 0}, Llf[4] = {0, 0, 0, 0}, Llgad[4] = {0, 0, 0, 0}, Llf2[4] = {0, 0, 0, 0};
        for (int d : days) {
            const int q = smp::PeriodOf(d);
            Lan[q] += La[d]; Lans[q] += Las[d]; Lall[q] += dl[d].L6;
            Llf[q] += La[d] * dl[d].LF32 * smp::kPT3fac;    // L x (LF3/LF2 x 64/63)
            Llf2[q] += dl[d].L6 * dl[d].LF32 * smp::kPT3fac;
            Llgad[q] += dl[d].Lstd_auth;
        }
        double T[6] = {0, 0, 0, 0, 0, 0};
        for (int q = 0; q < 4; ++q) { T[0] += Lan[q]; T[1] += Lans[q]; T[2] += Lall[q]; T[3] += Llf[q]; T[4] += Llgad[q]; T[5] += Llf2[q]; }
        FILE* o = fopen(Form("%s/smp_lumi.txt", smp::kTabOut), "w");
        fprintf(o, "# smp_lumi.txt -- written by plots/exp_report/smp_report.C (%s)\n", TDatime().AsString());
        fprintf(o, "# START-free (LGAD-free) luminosity from the HADES elastic-scattering analysis (talk LGADfree_luminosity_talk_draft_v2,\n");
        fprintf(o, "#   files of revision N, total L6 = 6080.5 nb^-1; copies in %s)\n", smp::kLumiDir);
        fprintf(o, "# L6_SF = six sectors (recommended); errors: stat, sectors (+ sector-2 correction; NOT independent between periods; a part of the\n"
                   "#   HADES syst), HADES syst (total, without sigma_int and the Argonne normalisation)\n");
        fprintf(o, "# L_an = sum of L6_SF over the runs present in the data (EpEm_ID rows); L_an(scaled) = with the two partly merged runs scaled\n");
        fprintf(o, "# f_LF = PT3/PT2 live fraction LF3/LF2 x 64/63 (L-weighted over the days); L x f_LF = luminosity if the live-fraction ratio applies\n");
        fprintf(o, "# L_LGAD = published chain with the LGAD START condition (Lstd_auth); L5 = five sectors x 6/5 (published convention)\n");
        fprintf(o, "%-5s %-8s | %9s %7s %7s %8s | %9s %9s | %9s %6s | %9s %9s\n", "per", "days", "L6_SF", "stat", "sect", "HADES", "L_an", "L_an(sc)",
                "L_an*fLF", "f_LF", "L_LGAD", "L5_SF");
        for (int q = 0; q < 4; ++q) {
            const smp::PerLumi& x = pl[q];
            fprintf(o, "%-5s %-8s | %9.3f %7.3f %7.3f %8.3f | %9.3f %9.3f | %9.3f %6.4f | %9.3f %9.3f\n", x.name.c_str(), smp::PeriodDays(q), x.L6, x.L6stat,
                    x.L6unc, x.L6hades, Lan[q], Lans[q], Llf[q], Llf[q] / Lan[q], x.Lstd_auth, x.L5);
            if (std::fabs(x.L6 - Lall[q]) > 0.01) fprintf(o, "#  WARNING: per-day sum %.3f differs from the per-period file\n", Lall[q]);
        }
        const smp::PerLumi& a = pl[4];
        fprintf(o, "%-5s %-8s | %9.3f %7.3f %7.3f %8.3f | %9.3f %9.3f | %9.3f %6.4f | %9.3f %9.3f\n", a.name.c_str(), "29 days", a.L6, a.L6stat, a.L6unc,
                a.L6hades, T[0], T[1], T[3], T[3] / T[0], a.Lstd_auth, a.L5);
        const double relH = smp::kLtotHades / smp::kLtot, relS = smp::kLtotSigInt / smp::kLtot, relN = smp::kLtotNorm / smp::kLtot;
        fprintf(o, "\n# total L6 = %.1f +- %.1f (stat) +- %.1f (HADES syst, %.2f %%) +- %.1f (sigma_int fit, %.2f %%) +- %.1f (Argonne norm, %.1f %%) nb^-1;"
                   " total %.2f %%\n", smp::kLtot, smp::kLtotStat, smp::kLtotHades, 100 * relH, smp::kLtotSigInt, 100 * relS, smp::kLtotNorm, 100 * relN,
                100 * std::sqrt(relH * relH + relS * relS + relN * relN));
        fprintf(o, "# analysed runs: %.3f nb^-1 (%.4f %% below all runs; %d of %d runs; the %d missing runs have no reconstructed particles in the\n"
                   "#   lepton production, %d of them with L = 0), scaled for the two partly merged runs %.3f\n",
                T[0], 100 * (1 - T[0] / a.L6), (int)std::count_if(rp.begin(), rp.end(), [](const std::pair<const int, smp::RunRows>& p) { return p.second.epem > 0; }),
                (int)rl.size(), (int)rl.size() - (int)std::count_if(rp.begin(), rp.end(), [](const std::pair<const int, smp::RunRows>& p) { return p.second.epem > 0; }),
                (int)std::count_if(rl.begin(), rl.end(), [&](const std::pair<const int, smp::RunLumi>& p) {
                    auto it = rp.find(p.first); return (it == rp.end() || it->second.epem == 0) && p.second.L6 == 0; }), T[1]);
        fprintf(o, "# PT3 rules: visible yield = N x 64/63 (trigger word 8192 = 63/64 of the PT3 events);\n");
        fprintf(o, "#   if the PT3/PT2 live-fraction ratio applies: L x f_LF = %.1f nb^-1 (f_LF = %.4f, per day %.4f-%.4f), i.e. L lower by %.1f %% (one-sided)\n",
                T[3], T[3] / T[0], [&]() { double m = 9; for (int d : days) m = std::min(m, dl[d].LF32 * smp::kPT3fac); return m; }(),
                [&]() { double m = 0; for (int d : days) m = std::max(m, dl[d].LF32 * smp::kPT3fac); return m; }(), 100 * (1 - T[3] / T[0]));
        fprintf(o, "# context: LGAD START condition (publication chain, six sectors) %.1f nb^-1 (= %.3f x L6_SF); after a START-finder repair (not final) %.0f nb^-1\n",
                smp::kLlgad, smp::kLlgad / smp::kLtot, smp::kLrepaired);
        fprintf(o, "# sector 2, days 034-044: hardware trigger deficit NOT contained in L (S2fac6 per day %.3f-%.3f is the correction of the PT2 elastic count)\n",
                [&]() { double m = 9; for (int d : days) if (smp::Sec2Deficit(d)) m = std::min(m, dl[d].S2fac6); return m; }(),
                [&]() { double m = 0; for (int d : days) if (smp::Sec2Deficit(d)) m = std::max(m, dl[d].S2fac6); return m; }());
        // no_hf runs
        fprintf(o, "\n# runs with no_hf = 1 (no elastic event in the window, L = 0) that are present in the data:\n");
        fprintf(o, "#   %-14s %8s %4s | %8s | %8s %8s %8s %8s %8s %8s\n", "dst_stem", "run", "day", "EpEm_ID", "us PAT", "us FAT", "pp PAT", "pp FAT", "mm PAT", "mm FAT");
        long long sumE = 0, sumU = 0, sumU3 = 0;
        for (const auto& p : rl) {
            if (!p.second.noHF) continue;
            auto it = rp.find(p.first);
            if (it == rp.end() || it->second.epem == 0) continue;
            const smp::RunRows& r = it->second;
            fprintf(o, "#   %-14s %8d %4d | %8lld | %8lld %8lld %8lld %8lld %8lld %8lld\n", p.second.stem.c_str(), p.first, p.second.day, r.epem, r.us0, r.us3,
                    r.pp0, r.pp3, r.mm0, r.mm3);
            sumE += r.epem; sumU += r.us0; sumU3 += r.us3;
        }
        long long totE = 0, totU3 = 0; for (const auto& p : rp) { totE += p.second.epem; totU3 += p.second.us3; }
        fprintf(o, "#   sum: EpEm_ID %lld of %lld rows (%.4f %%); FAT e+e- rows %lld of %lld (%.4f %%); after pre e+e- %.0f, e+e+ %.0f, e-e- %.0f rows\n",
                sumE, totE, 100. * sumE / totE, sumU3, totU3, 100. * sumU3 / totU3, nohf[1][0], nohf[1][1], nohf[1][2]);
        fprintf(o, "#   effect of removing them on the eta peak (all days): see smp_peak_periods.txt (pre S 4670.6 -> 4670.8): not material\n");
        fclose(o);
        P("lumi: L_an %.3f (scaled %.3f), L6 all %.3f, L x f_LF %.3f (f %.4f), L6xLF (all runs) %.3f, LGAD per-day sum %.3f\n", T[0], T[1], T[2], T[3], T[3] / T[0], T[5], T[4]);

        FILE* x = fopen(Form("%s/smp_lumi.tex", smp::kTabOut), "w");
        fprintf(x, "%% smp_lumi.tex -- written by plots/exp_report/smp_report.C (numbers: tables_exp/smp_lumi.txt)\n");
        fprintf(x, "\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{4pt}\n");
        fprintf(x, "\\caption{Integrated luminosity per period in nb$^{-1}$. $L_\\mathrm{SF}$: START-free luminosity of the HADES elastic-scattering "
                   "analysis, six sectors, all runs, with the statistical error and the sector uncertainty (the latter is a part of the 5.56\\,\\%% HADES "
                   "systematic and is not independent between periods). "
                   "$L_\\mathrm{an}$: sum over the runs present in the data (17\\,557 of 17\\,615). $f_\\mathrm{LF}$: PT3/PT2 live-fraction ratio "
                   "$\\mathrm{LF3}/\\mathrm{LF2}\\times 64/63$, $L$-weighted; if it applies, $L$ is lower by up to $1-f_\\mathrm{LF}$ (one-sided). "
                   "$L_\\mathrm{LGAD}$: published chain with the LGAD START condition. Total $L_\\mathrm{SF}$ uncertainty: "
                   "$\\pm5.56\\,\\%%$ (HADES) $\\pm1.03\\,\\%%$ ($\\sigma_\\mathrm{int}$ fit) $\\pm4.0\\,\\%%$ (normalisation), %.1f\\,\\%% in total. "
                   "The errors of the last line are taken from the luminosity file: they include correlations between the periods (e.g.\\ a 100\\,\\%%-correlated "
                   "simulation term) and are not the quadratic or linear sums of the period values.}\n",
                100 * std::sqrt(relH * relH + relS * relS + relN * relN));
        fprintf(x, "\\label{tab:smp_lumi}\n\\begin{tabular}{ll rr rr r}\n\\toprule\n");
        fprintf(x, "Period & Days & $L_\\mathrm{SF}$ & $L_\\mathrm{an}$ & $f_\\mathrm{LF}$ & $L_\\mathrm{an} f_\\mathrm{LF}$ & $L_\\mathrm{LGAD}$ \\\\\n\\midrule\n");
        for (int q = 0; q < 4; ++q)
            fprintf(x, "%s & %s & $%.1f\\pm%.1f\\pm%.1f$ & %.1f & %.4f & %.1f & %.1f \\\\\n", pl[q].name.c_str(), smp::PeriodDays(q), pl[q].L6, pl[q].L6stat,
                    pl[q].L6unc, Lan[q], Llf[q] / Lan[q], Llf[q], pl[q].Lstd_auth);
        fprintf(x, "\\midrule\nall & 29 days & $%.1f\\pm%.1f\\pm%.1f$ & %.1f & %.4f & %.1f & %.1f \\\\\n", a.L6, a.L6stat, a.L6unc, T[0], T[3] / T[0], T[3], a.Lstd_auth);
        fprintf(x, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        fclose(x);
    }

    // ---------------- 3. cut flow ----------------
    {
        FILE* o = fopen(Form("%s/smp_cutflow.txt", smp::kTabOut), "w");
        fprintf(o, "# smp_cutflow.txt -- written by plots/exp_report/smp_report.C (%s)\n", TDatime().AsString());
        fprintf(o, "# full cut flow of the three hypotheses; rows = stored particle combinations, events = event groups (adjacent rows with equal\n");
        fprintf(o, "# (file, run, event)); steps 0-3 from the PAT files (smp_pat_counts.C), steps from 'FAT rows' on from the compact file (smp_flow.C).\n");
        fprintf(o, "# Steps are cumulative down to 'dz_pi < 40'; the named selections below the line are each defined from pre.\n");
        fprintf(o, "# frac = rows kept relative to the previous cumulative step (named selections: relative to pre).\n");
        fprintf(o, "# isBest == 1 rows = events of the isBest == 1 variant. OA = OA(e+e-) > 2 deg variant (not a fixed cut).\n\n");
        long long pr[3][4], pe[3][4], pb[3][4];
        for (int k = 0; k < 3; ++k) for (int s = 0; s < 4; ++s) { pr[k][s] = pe[k][s] = pb[k][s] = 0; for (int d : days) { pr[k][s] += pat[d].rows[k][s]; pe[k][s] += pat[d].ev[k][s]; pb[k][s] += pat[d].best[k][s]; } }
        // consistency: PAT step 3 = FAT rows
        for (int k = 0; k < 3; ++k)
            fprintf(o, "# check %s: PAT rows after PT3 %lld, FAT rows %.0f; events %lld vs %.0f; isBest==1 %lld vs %.0f\n", hn[k], pr[k][3], C(0, 0, k).rows,
                    pe[k][3], C(0, 0, k).ev, pb[k][3], C(1, 0, k).rows);
        struct L { std::string name; double r[3], e[3], b[3], ro[3]; bool chain; };
        std::vector<L> lines;
        const char* patName[4] = {"PAT rows", "isBest >= 0 (FAT isBestAll)", "vz > -500 mm", "PT3 (trigbit == 8192) = FAT rows"};
        for (int s = 0; s < 4; ++s) {
            L l; l.name = patName[s]; l.chain = true;
            for (int k = 0; k < 3; ++k) { l.r[k] = pr[k][s]; l.e[k] = pe[k][s]; l.b[k] = pb[k][s]; l.ro[k] = -1; }
            lines.push_back(l);
        }
        const char* fatName[13] = {"", "cor_ok == 1", "vz > -160 mm (= pre = BUMP)", "both pions in the vertex fit", "dz_pi < 40 mm",
                                   "DIRECT = pre + dz_pi < 40 mm", "PEAK = pre + dz_pi < 40 + MM4 < 2.85", "peak nominal", "peak nominal + vertex",
                                   "bump nominal", "bump nominal + vertex", "X17 region (sim window)", "X17 region (data window)"};
        for (int i = 1; i < 13; ++i) {
            L l; l.name = fatName[i]; l.chain = (i <= 4);
            for (int k = 0; k < 3; ++k) { l.r[k] = C(0, i, k).rows; l.e[k] = C(0, i, k).ev; l.b[k] = C(1, i, k).rows; l.ro[k] = C(2, i, k).rows; }
            lines.push_back(l);
        }
        fprintf(o, "%-38s |", "step");
        for (int k = 0; k < 3; ++k) fprintf(o, " %9s %6s %9s %9s %9s |", Form("%s rows", hn[k]), "frac", "events", "isBest==1", "OA>2 rows");
        fprintf(o, "\n");
        const int iPre = 4 + 1;    // index of pre in 'lines'
        for (size_t j = 0; j < lines.size(); ++j) {
            const L& l = lines[j];
            if (j == 8) fprintf(o, "%s\n", std::string(38 + 3 * 52, '-').c_str());
            fprintf(o, "%-38s |", l.name.c_str());
            for (int k = 0; k < 3; ++k) {
                double ref = 0;
                if (j > 0 && l.chain) ref = lines[j - 1].r[k];
                else if (!l.chain) ref = lines[iPre].r[k];
                char ro[16]; if (l.ro[k] >= 0) snprintf(ro, 16, "%.0f", l.ro[k]); else snprintf(ro, 16, "-");
                fprintf(o, " %9.0f %6.4f %9.0f %9.0f %9s |", l.r[k], ref > 0 ? l.r[k] / ref : 1., l.e[k], l.b[k], ro);
            }
            fprintf(o, "\n");
        }
        // vertex cut detail
        fprintf(o, "\n# vertex cut vz > -500 mm at the PAT level (rows with isBest >= 0):\n");
        for (int k = 0; k < 3; ++k) {
            long long nv = 0, nf = 0; for (int d : days) { nv += pat[d].vzNoVtx[k]; nf += pat[d].vzFail[k]; }
            double mn = 1, mx = 0;
            for (int d : days) { const double f = 1.0 - (double)pat[d].rows[k][2] / pat[d].rows[k][1]; mn = std::min(mn, f); mx = std::max(mx, f); }
            fprintf(o, "#   %s: removed %lld of %lld rows (%.1f %%; per day %.1f-%.1f %%); of the removed rows %.2f %% have vz = -1000 (no vertex reconstructed)\n",
                    hn[k], nf, pr[k][1], 100. * nf / pr[k][1], 100 * mn, 100 * mx, 100. * nv / nf);
            P("vertex cut %s: removed %.2f %% (days %.1f-%.1f %%), no-vertex share %.2f %%\n", hn[k], 100. * nf / pr[k][1], 100 * mn, 100 * mx, 100. * nv / nf);
        }
        fprintf(o, "# trigger words after the vertex cut (rows): ");
        for (int k = 0; k < 3; ++k) {
            std::map<int, long long> t; for (int d : days) for (auto& q : pat[d].trigAfter[k]) t[q.first] += q.second;
            fprintf(o, "%s:", hn[k]);
            for (auto& q : t) fprintf(o, " %d %lld", q.first, q.second);
            fprintf(o, " ; ");
        }
        fprintf(o, "\n# (FAT keeps trigbit == 8192 exactly; rows with 12288 = 8192 + 4096 are dropped)\n");
        long long nBN = 0, nG1 = 0; for (int d : days) for (int k = 0; k < 3; ++k) { nBN += pat[d].isBestNeg[k]; nG1 += pat[d].grpBestNot1[k]; }
        fprintf(o, "# PAT checks: rows with isBest < 0: %lld; event groups without exactly one isBest == 1 row: %lld\n", nBN, nG1);

        // windows (e+e-)
        fprintf(o, "\n# e+e- rows in the mass windows and CB/N (CB = 2 sqrt(N++ N--), 10 MeV bins interpolated, k = 1):\n");
        fprintf(o, "#   pk = data peak 0.53364-0.55416, bu = data bump 0.4550-0.5225, pkS/buS = frozen simulation windows 0.5333-0.5587 / 0.4576-0.5255 GeV/c2\n");
        for (int v = 0; v < 4; ++v) {
            const char* vn[4] = {"isBest >= 0", "isBest == 1", "isBest >= 0, OA > 2", "isBest == 1, OA > 2"};
            fprintf(o, "\n# variant %s\n%-38s | %8s %6s | %8s %6s | %8s %6s | %8s %6s\n", vn[v], "step", "N pk", "CB/N", "N bu", "CB/N", "N pkS", "CB/N", "N buS", "CB/N");
            for (int i = 1; i < 13; ++i) {
                const CF& c = C(v, i, 0);
                fprintf(o, "%-38s |", fatName[i]);
                for (int w = 0; w < 4; ++w) {
                    // X17 rows: no bump-window entries; the data-window X17 row has no meaning in the wider sim window
                    const bool undef = (i >= 11 && (w == 1 || w == 3)) || (i == 12 && w == 2);
                    if (undef) fprintf(o, " %8s %6s |", "-", "-");
                    else fprintf(o, " %8.0f %6.4f |", c.w[w], c.w[w] > 0 ? c.cb[w] / c.w[w] : 0.);
                }
                fprintf(o, "\n");
            }
        }
        fclose(o);

        // LaTeX: main table
        FILE* x = fopen(Form("%s/smp_cutflow.tex", smp::kTabOut), "w");
        fprintf(x, "%% smp_cutflow.tex -- written by plots/exp_report/smp_report.C (numbers: tables_exp/smp_cutflow.txt)\n");
        fprintf(x, "\\begin{table}[H]\n\\centering\\scriptsize\n\\setlength{\\tabcolsep}{3pt}\n");
        fprintf(x, "\\caption{Cut flow of the three hypotheses, all 29 days. Rows = stored particle combinations (all combinations, isBest $\\geq 0$), "
                   "events = event groups; the column isBest $=1$ counts the best combinations (= events of that variant). The steps above the line are "
                   "cumulative; each named selection below it is defined from pre: DIRECT = pre + $\\Delta z_\\pi<40$~mm, PEAK = DIRECT + $MM_4<2.85$~GeV/$c^2$ "
                   "(frozen selections of the simulation study; pre is also the frozen BUMP selection), peak nominal = pre + $\\alpha_\\mathrm{lab}<55^\\circ$ + "
                   "$M_{\\pi\\pi}>0.30$~GeV/$c^2$, bump nominal = pre + $\\alpha_\\mathrm{lab}<40^\\circ$ + $0.30<M_{\\pi\\pi}<0.41$~GeV/$c^2$, "
                   "vertex = both pions in the event-vertex fit, $\\Delta z_\\pi=\\max(|z_{\\pi^+}-v_z|,|z_{\\pi^-}-v_z|)$. The first four lines are counted in the PAT files, the others in the "
                   "FAT output. Corrected momenta (cor\\_ok $=1$) from that step on. The frozen windows of the simulation study are used for the X17 region "
                   "(0.5333--0.5587~GeV/$c^2$, $15.5\\leq M_{ee}<24.5$~MeV/$c^2$); with the data peak window 0.53364--0.55416~GeV/$c^2$ the numbers are "
                   "%.0f, %.0f and %.0f rows. Not normalised to luminosity.}\n", C(0, 12, 0).rows, C(0, 12, 1).rows, C(0, 12, 2).rows);
        fprintf(x, "\\label{tab:smp_cutflow}\n\\begin{tabular}{l rrr rrr rrr}\n\\toprule\n");
        fprintf(x, " & \\multicolumn{3}{c}{$\\pi^+\\pi^-e^+e^-$} & \\multicolumn{3}{c}{$\\pi^+\\pi^-e^+e^+$} & \\multicolumn{3}{c}{$\\pi^+\\pi^-e^-e^-$} \\\\\n");
        fprintf(x, "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\n");
        fprintf(x, "Step & rows & events & isBest$=1$ & rows & events & isBest$=1$ & rows & events & isBest$=1$ \\\\\n\\midrule\n");
        const char* texName[17] = {"PAT rows", "isBest $\\geq 0$", "$v_z>-500$~mm", "PT3 = FAT rows", "cor\\_ok $=1$",
                                   "$v_z>-160$~mm (pre)", "+ pions in vertex fit", "+ $\\Delta z_\\pi<40$~mm",
                                   "DIRECT", "PEAK", "peak nominal", "peak nominal + vertex",
                                   "bump nominal", "bump nominal + vertex", "X17 region", "", ""};
        for (size_t j = 0; j + 1 < lines.size(); ++j) {   // last line (X17 data window) in the caption
            if (j == 8) fprintf(x, "\\midrule\n");
            fprintf(x, "%s", texName[j]);
            for (int k = 0; k < 3; ++k) fprintf(x, " & %.0f & %.0f & %.0f", lines[j].r[k], lines[j].e[k], lines[j].b[k]);
            fprintf(x, " \\\\\n");
        }
        fprintf(x, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        fclose(x);

        // LaTeX: windows table (e+e-), without and with OA > 2
        FILE* y = fopen(Form("%s/smp_cutflow_windows.tex", smp::kTabOut), "w");
        fprintf(y, "%% smp_cutflow_windows.tex -- written by plots/exp_report/smp_report.C (numbers: tables_exp/smp_cutflow.txt)\n");
        fprintf(y, "\\begin{table}[H]\n\\centering\\footnotesize\n\\setlength{\\tabcolsep}{3.5pt}\n");
        fprintf(y, "\\caption{$e^+e^-$ rows in the data peak window (0.53364--0.55416~GeV/$c^2$) and the data bump window (0.4550--0.5225~GeV/$c^2$) "
                   "of $M(\\pi^+\\pi^-e^+e^-)$ and the combinatorial-background share CB/$N$ (CB $=2\\sqrt{N_{++}N_{--}}$ in 10~MeV bins, "
                   "interpolated, $k=1$), for all combinations (isBest $\\geq 0$), without and with OA$(e^+e^-)>2^\\circ$. "
                   "Steps as in Table~\\ref{tab:smp_cutflow}.}\n");
        fprintf(y, "\\label{tab:smp_cutflow_windows}\n\\begin{tabular}{l rr rr | rr rr}\n\\toprule\n");
        fprintf(y, " & \\multicolumn{4}{c|}{all OA} & \\multicolumn{4}{c}{OA $>2^\\circ$} \\\\\n");
        fprintf(y, "Step & $N_\\mathrm{peak}$ & CB/$N$ & $N_\\mathrm{bump}$ & CB/$N$ & $N_\\mathrm{peak}$ & CB/$N$ & $N_\\mathrm{bump}$ & CB/$N$ \\\\\n\\midrule\n");
        for (int i = 1; i < 12; ++i) {
            if (i == 5) fprintf(y, "\\midrule\n");
            fprintf(y, "%s", texName[i + 3]);
            for (int v : {0, 2}) {
                const CF& c = C(v, i, 0);
                if (i == 11) fprintf(y, " & %.0f & %.3f & -- & --", c.w[0], c.w[0] > 0 ? c.cb[0] / c.w[0] : 0.);   // X17: no bump-window rows
                else fprintf(y, " & %.0f & %.3f & %.0f & %.3f", c.w[0], c.w[0] > 0 ? c.cb[0] / c.w[0] : 0., c.w[1], c.w[1] > 0 ? c.cb[1] / c.w[1] : 0.);
            }
            fprintf(y, " \\\\\n");
        }
        fprintf(y, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        fclose(y);

        // ---------------- figure: cut flow ----------------
        TCanvas* c = new TCanvas("cCF", "", 1000, 640);
        c->SetLeftMargin(0.11); c->SetBottomMargin(0.36); c->SetRightMargin(0.03); c->SetTopMargin(0.05); c->SetLogy();
        const int nl = (int)lines.size() - 1;
        TH1D* f0 = new TH1D("frCF", "", nl, 0, nl);
        const char* rootName[16] = {"PAT counts", "isBest #geq 0", "v_{z} > -500 mm", "PT3 = FAT counts", "cor_ok", "v_{z} > -160 mm (pre)",
                                    "+ pions in vertex fit", "+ #Deltaz_{#pi} < 40 mm", "DIRECT", "PEAK", "peak nominal", "peak nom. + vertex",
                                    "bump nominal", "bump nom. + vertex", "X17 region", ""};
        for (int j = 0; j < nl; ++j) f0->GetXaxis()->SetBinLabel(j + 1, rootName[j]);
        f0->GetXaxis()->LabelsOption("v");
        f0->GetXaxis()->SetLabelSize(0.045);
        f0->SetMinimum(1e-5); f0->SetMaximum(3);
        f0->GetYaxis()->SetTitle("fraction of the PAT counts");
        f0->GetYaxis()->SetTitleOffset(0.95);
        f0->Draw("axis");
        TLine* sep = new TLine(8, 1e-5, 8, 3); sep->SetLineStyle(2); sep->SetLineColor(kGray + 1); sep->Draw();
        const int col[3] = {kBlack, kGreen + 2, kOrange + 7};   // [review fix] like-sign key of the report: e+e+ green, e-e- orange
        const int mk[3] = {20, 21, 22}, mko[3] = {24, 25, 26};
        TLegend* lg = new TLegend(0.13, 0.38, 0.56, 0.60);
        lg->SetNColumns(2); lg->SetTextSize(0.036); lg->SetFillStyle(0);
        for (int k = 0; k < 3; ++k) {
            TGraph* g = new TGraph(), *gb = new TGraph();
            for (int j = 0; j < nl; ++j) {
                g->SetPoint(j, j + 0.5 - 0.12 + 0.12 * k, lines[j].r[k] / lines[0].r[k]);
                gb->SetPoint(j, j + 0.5 - 0.12 + 0.12 * k, lines[j].b[k] / lines[0].r[k]);
            }
            g->SetMarkerStyle(mk[k]); g->SetMarkerColor(col[k]); g->SetLineColor(col[k]); g->SetMarkerSize(1.3);
            gb->SetMarkerStyle(mko[k]); gb->SetMarkerColor(col[k]); gb->SetLineColor(col[k]); gb->SetMarkerSize(1.3);
            g->Draw("P same"); gb->Draw("P same");
            lg->AddEntry(g, Form("#pi^{+}#pi^{-}%s, counts", k == 0 ? "e^{+}e^{-}" : (k == 1 ? "e^{+}e^{+}" : "e^{-}e^{-}")), "p");
            lg->AddEntry(gb, "isBest = 1", "p");
        }
        lg->Draw();
        TLatex t; t.SetTextSize(0.034); t.SetTextColor(kGray + 2);
        t.DrawLatex(8.2, 1.5e-5, "named selections (each from pre)");
        smp::Save(c, std::string(smp::kFigOut) + "/smp_cutflow");
    }

    // ---------------- 4. rates vs the live fraction ----------------
    {
        P("\n# per-day rates R (pre) against the PT3/PT2 live-fraction ratio f = LF3/LF2 x 64/63 and the PT2 live fraction LF2\n");
        for (int k = 0; k < 3; ++k) {
            std::vector<double> R, Rf, f, lf2;
            for (int d : days) { R.push_back(dr[d][k].R); f.push_back(dl[d].LF32 * smp::kPT3fac); lf2.push_back(dl[d].LF2); Rf.push_back(dr[d][k].R / f.back()); }
            auto rms = [](const std::vector<double>& v) { double m = 0, s = 0; for (double x : v) m += x; m /= v.size(); for (double x : v) s += (x / m - 1) * (x / m - 1); return std::sqrt(s / v.size()); };
            auto corr = [](const std::vector<double>& a, const std::vector<double>& b) {
                const int n = (int)a.size(); double ma = 0, mb = 0; for (int i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; } ma /= n; mb /= n;
                double sab = 0, saa = 0, sbb = 0;
                for (int i = 0; i < n; ++i) { sab += (a[i] - ma) * (b[i] - mb); saa += (a[i] - ma) * (a[i] - ma); sbb += (b[i] - mb) * (b[i] - mb); }
                return sab / std::sqrt(saa * sbb); };
            const double c1 = corr(R, f), c2 = corr(R, lf2);
            P("%s: rms(R/<R>) = %.4f ; rms after dividing by f: %.4f ; corr(R, f) = %+.3f ; corr(R, LF2) = %+.3f\n", hn[k], rms(R), rms(Rf), c1, c2);
        }
    }

    // ---------------- 5. figure: vertex cut ----------------
    {
        TH1D* hv[3] = {nullptr, nullptr, nullptr};
        const char* tag[3] = {"us", "pp", "mm"};
        for (int d : days) {
            TFile* f = TFile::Open(Form("%s/pat/pat_%03d.root", smp::kWork, d), "READ");
            for (int k = 0; k < 3; ++k) {
                TH1D* h = (TH1D*)f->Get(Form("vz_%s", tag[k]));
                if (!hv[k]) { hv[k] = (TH1D*)h->Clone(Form("hv%d", k)); hv[k]->SetDirectory(nullptr); }
                else hv[k]->Add(h);
            }
            f->Close(); delete f;
        }
        TCanvas* c = new TCanvas("cV", "", 1500, 520);
        c->Divide(3, 1, 0.002, 0.002);
        const int col[3] = {kBlack, kGreen + 2, kOrange + 7};   // [review fix] like-sign key of the report: e+e+ green, e-e- orange
        const int mk[3] = {20, 21, 22};
        const char* hyp[3] = {"#pi^{+}#pi^{-}e^{+}e^{-}", "#pi^{+}#pi^{-}e^{+}e^{+}", "#pi^{+}#pi^{-}e^{-}e^{-}"};
        auto axes = [](TH1* h) {
            h->GetXaxis()->SetLabelSize(0.056); h->GetXaxis()->SetTitleSize(0.058); h->GetXaxis()->SetTitleOffset(1.05);
            h->GetYaxis()->SetLabelSize(0.056); h->GetYaxis()->SetTitleSize(0.058); h->GetYaxis()->SetTitleOffset(1.30);
            h->GetXaxis()->SetNdivisions(505);
        };
        // (a) full range, log y, normalised to all rows
        TPad* p1 = (TPad*)c->cd(1);
        p1->SetLogy(); p1->SetLeftMargin(0.17); p1->SetRightMargin(0.03);
        TLegend* lg = new TLegend(0.20, 0.68, 0.97, 0.90);
        lg->SetTextSize(0.048); lg->SetFillStyle(0);
        TH1D* hz[3]; TH1D* h1[3];
        for (int k = 0; k < 3; ++k) {
            hz[k] = (TH1D*)hv[k]->Clone(Form("hz%d", k)); hz[k]->SetDirectory(nullptr);   // 1 mm bins, for the zoom
            h1[k] = (TH1D*)hv[k]->Clone(Form("h1%d", k)); h1[k]->SetDirectory(nullptr);   // 1 mm bins, for the numbers
            TH1D* h = hv[k];
            h->Rebin(10);
            const double tot = h->Integral(0, h->GetNbinsX() + 1);
            h->Scale(1.0 / tot);
            h->SetLineColor(col[k]); h->SetLineWidth(2);
            h->GetXaxis()->SetRangeUser(-1050, 300);
            h->SetMinimum(1e-6); h->SetMaximum(30);
            h->GetXaxis()->SetTitle("event vertex v_{z} [mm]");
            h->GetYaxis()->SetTitle("fraction of the PAT counts per 10 mm");
            axes(h);
            h->Draw(k == 0 ? "hist" : "hist same");
            const double fl = h->Integral(h->FindBin(-1000.0), h->FindBin(-1000.0));
            lg->AddEntry(h, Form("%s: %.1f %% at v_{z} = -1000 mm", hyp[k], 100 * fl), "l");
        }
        lg->Draw();
        { TLatex tl; tl.SetNDC(); tl.SetTextSize(0.058); tl.DrawLatex(0.18, 0.93, "(a) full range"); }
        TLine* l1 = new TLine(-500, 1e-6, -500, 0.05); l1->SetLineStyle(2); l1->SetLineColor(kBlue + 1); l1->SetLineWidth(2); l1->Draw();
        TLine* l2 = new TLine(-160, 1e-6, -160, 0.05); l2->SetLineStyle(3); l2->SetLineColor(kBlue + 1); l2->SetLineWidth(2); l2->Draw();
        TLatex t; t.SetTextSize(0.042); t.SetTextColor(kBlue + 1); t.SetTextAlign(21);
        t.DrawLatex(-500, 0.08, "-500"); t.DrawLatex(-160, 0.08, "-160");
        t.SetTextColor(kGray + 2); t.SetTextSize(0.040); t.SetTextAlign(12);
        t.DrawLatex(-975, 0.25, "no vertex");
        // (b) zoom on the target region, normalised to the rows with vz > -500
        TPad* p2 = (TPad*)c->cd(2);
        p2->SetLeftMargin(0.17); p2->SetRightMargin(0.03);
        double ymx = 0;
        for (int k = 0; k < 3; ++k) {
            TH1D* h = hz[k];
            h->Rebin(4);
            const double tot = h->Integral(h->FindBin(-499.9), h->GetNbinsX() + 1);
            h->Scale(1.0 / tot);
            h->GetXaxis()->SetRangeUser(-260, 60);
            ymx = std::max(ymx, h->GetMaximum());
        }
        // vz regions (rows with vz > -500 mm; 1 mm histograms, bin edges on integer mm)
        {
            const double ed[8] = {-500, -200, -180, -170, -160, -150, -100, 1e9};
            P("\n# vz regions of the PAT rows with vz > -500 mm (all days, isBest >= 0): fraction per region [mm]\n#   %-26s", "");
            for (int r = 0; r < 7; ++r) P(" %13s", r < 6 ? Form("(%.0f,%.0f]", ed[r], ed[r + 1]) : "> -100");
            P("\n");
            for (int k = 0; k < 3; ++k) {
                TH1D* h = h1[k];
                const double tot = h->Integral(h->FindBin(-499.5), h->GetNbinsX() + 1);
                P("#   %-26s", hyp[k]);
                for (int r = 0; r < 7; ++r) {
                    const int b1 = h->FindBin(ed[r] + 0.5), b2 = r < 6 ? h->FindBin(ed[r + 1] - 0.5) : h->GetNbinsX() + 1;
                    P(" %13.4f", h->Integral(b1, b2) / tot);
                }
                P("\n");
                int bm = 1; double vm = 0;
                for (int b = h->FindBin(-300.); b <= h->FindBin(100.); ++b) {
                    const double v = h->GetBinContent(b - 1) + h->GetBinContent(b) + h->GetBinContent(b + 1);
                    if (v > vm) { vm = v; bm = b; }
                }
                P("#     mode at vz = %.1f mm (3 mm running sum, 1 mm bins)\n", h->GetBinCenter(bm));
            }
        }
        TLegend* lgz = new TLegend(0.20, 0.72, 0.97, 0.92);
        lgz->SetTextSize(0.048); lgz->SetFillStyle(0);
        for (int k = 0; k < 3; ++k) {
            TH1D* h = hz[k];
            h->SetLineColor(col[k]); h->SetLineWidth(2);
            h->SetMinimum(0); h->SetMaximum(1.45 * ymx);
            h->GetXaxis()->SetTitle("event vertex v_{z} [mm]");
            h->GetYaxis()->SetTitle("fraction of counts (v_{z} > -500) per 4 mm");
            axes(h); h->GetYaxis()->SetTitleSize(0.054);
            h->Draw(k == 0 ? "hist" : "hist same");
            const double fk = h->Integral(h->FindBin(-159.99), h->GetNbinsX() + 1);
            lgz->AddEntry(h, Form("%s: %.3f with v_{z} > -160 mm", hyp[k], fk), "l");
        }
        lgz->Draw();
        { TLatex tl; tl.SetNDC(); tl.SetTextSize(0.058); tl.DrawLatex(0.18, 0.93, "(b) target region"); }
        TLine* l3 = new TLine(-160, 0, -160, 1.0 * ymx); l3->SetLineStyle(3); l3->SetLineColor(kBlue + 1); l3->SetLineWidth(2); l3->Draw();
        // (c) fraction kept per day
        TPad* p3 = (TPad*)c->cd(3);
        p3->SetLeftMargin(0.17); p3->SetRightMargin(0.03);
        TH1D* f0 = new TH1D("frV", "", 34, 33.5, 67.5);
        f0->SetMinimum(0); f0->SetMaximum(1.0);
        f0->GetXaxis()->SetTitle("day of 2022");
        f0->GetYaxis()->SetTitle("fraction of counts kept by v_{z} > -500 mm");
        axes(f0); f0->GetYaxis()->SetTitleSize(0.054);
        f0->Draw("axis");
        TLegend* lg2 = new TLegend(0.20, 0.40, 0.90, 0.58);
        lg2->SetTextSize(0.048); lg2->SetFillStyle(0);
        for (int k = 0; k < 3; ++k) {
            TGraph* g = new TGraph();
            int n = 0;
            for (int d : days) g->SetPoint(n++, d, (double)pat[d].rows[k][2] / pat[d].rows[k][1]);
            g->SetMarkerStyle(mk[k]); g->SetMarkerColor(col[k]); g->SetLineColor(col[k]); g->SetMarkerSize(1.1);
            g->Draw("P same");
            long long a = 0, b = 0; for (int d : days) { a += pat[d].rows[k][2]; b += pat[d].rows[k][1]; }
            lg2->AddEntry(g, Form("%s: %.3f kept (all days)", hyp[k], (double)a / b), "p");
        }
        lg2->Draw();
        { TLatex tl; tl.SetNDC(); tl.SetTextSize(0.058); tl.DrawLatex(0.18, 0.93, "(c) kept by v_{z} > #minus500 mm per day"); }
        smp::Save(c, std::string(smp::kFigOut) + "/smp_vertex_cut");
    }

    fclose(gOut); gOut = nullptr;
}
