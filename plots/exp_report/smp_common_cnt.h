// smp_common_cnt.h -- COPY of plots/exp_report/smp_common.h made by figure agent A (integ4, 2026-11-15) with
// the figure strings 'rows' -> 'counts' (analyst's terminology) and the outputs of the copy redirected:
// figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/smp
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/mk_cnt.py; the physics code is unchanged.
// ============================================================================
// smp_common.h  --  shared constants and readers of the chapter "samples,
//                   luminosity and cut flow" (prefix smp) of the data report
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Header only, included by plots/exp_report/smp_*.C. It holds
//   - the paths of the inputs and outputs of the chapter,
//   - the 29 days of the beam time and their periods
//     I (034-046), IIA (047-048), IIB (049-053), III (059-067), and the
//     sector-2 trigger-deficit days 034-044,
//   - readers of the START-free luminosity tables (revision N, total
//     6080.5 nb^-1; copies of the files of the HADES elastic-scattering
//     analysis in <work>/lumi_inputs/):
//       sf_luminosity_perperiod.csv  per period
//       lumi_startfree_perday.csv    per day (L6, live fractions LF2, LF3, ...)
//       lumi_coverage_perrun.csv     per run, with the number of DSTs of the run
//                                    that are in the FAT input (n_dst_merged)
//   - the plot style and the PDF/PNG writer (vector PDF via EPS + ps2pdf -dEPSCrop).
// No run command: it is not a macro.
// ============================================================================
#ifndef SMP_COMMON_CNT_H
#define SMP_COMMON_CNT_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

namespace smp {

// ---- paths -------------------------------------------------------------------
const char* const kPatList  = "/hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list";
const char* const kPatDir   = "/hdd3/przygoda/HADES/files/exp/gen4";
const char* const kProd     = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530";
const char* const kWork     = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/smp";
const char* const kLumiDir  = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/smp/lumi_inputs";
const char* const kFigOut   = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp";
const char* const kTabOut   = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/smp/tab";
const char* const kWorkOut  = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/smp/work";      // written files of the copy (flow.txt, report.txt)
const char* const kCntScratch = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/smp/figs";
const char* const kCntScope = " smp_cutflow smp_vertex_cut smp_peak_periods smp_peak_spectra_nomvtx smp_peak_spectra_pre smp_rates_perday ";

// ---- PT3 normalisation ---------------------------------------------------------
// trigger word 8192 holds 63/64 of the PT3 events: visible yield = N x 64/63
const double kPT3fac  = 64.0 / 63.0;
// PT3-input live fraction relative to PT2 (L-weighted, incl. 64/63): 0.932
const double kLF32    = 0.932;

// ---- luminosity numbers quoted from the luminosity analysis (nb^-1) ------------
const double kLtot        = 6080.524;   // L6_SF, 17 615 runs, six sectors
const double kLtotStat    = 5.172;
const double kLtotHades   = 337.942;    // 5.56 %
const double kLtotSigInt  = 62.6;       // sigma_int (Argonne) fit, 1.03 %
const double kLtotNorm    = 243.2;      // Argonne normalisation 4.0 %
const double kLlgad       = 5587.53;    // published chain with the LGAD START condition (Lstd_auth)
const double kLrepaired   = 6722.0;     // after a START-finder repair (not final), quoted from the talk

// ---- days and periods ----------------------------------------------------------
inline std::vector<int> Days() {
    std::vector<int> d;
    for (int i = 34; i <= 53; ++i) d.push_back(i);
    for (int i = 59; i <= 67; ++i) d.push_back(i);
    return d;
}
// period index 0 = I, 1 = IIA, 2 = IIB, 3 = III
inline int PeriodOf(int day) {
    if (day >= 34 && day <= 46) return 0;
    if (day >= 47 && day <= 48) return 1;
    if (day >= 49 && day <= 53) return 2;
    if (day >= 59 && day <= 67) return 3;
    return -1;
}
const int kNPer = 4;
inline const char* PeriodName(int p) {
    static const char* n[4] = {"I", "IIA", "IIB", "III"};
    return (p >= 0 && p < 4) ? n[p] : "?";
}
inline const char* PeriodDays(int p) {
    static const char* n[4] = {"034-046", "047-048", "049-053", "059-067"};
    return (p >= 0 && p < 4) ? n[p] : "?";
}
// sector-2 hardware trigger deficit (not in L): days 034-044
inline bool Sec2Deficit(int day) { return day >= 34 && day <= 44; }

// ---- small CSV reader ------------------------------------------------------------
// returns the rows (split at commas, quotes removed) after the header line that
// starts with 'firstcol'; lines starting with '#' are skipped
struct Csv {
    std::vector<std::string> head;
    std::vector<std::vector<std::string>> rows;
    int col(const char* name) const {
        for (size_t i = 0; i < head.size(); ++i) if (head[i] == name) return (int)i;
        printf("smp::Csv: column %s missing\n", name);
        return -1;
    }
};
inline std::vector<std::string> SplitCsv(const std::string& line) {
    std::vector<std::string> out; std::string cur; bool q = false;
    for (char ch : line) {
        if (ch == '"') { q = !q; continue; }
        if (ch == ',' && !q) { out.push_back(cur); cur.clear(); continue; }
        if (ch == '\r' || ch == '\n') continue;
        cur += ch;
    }
    out.push_back(cur);
    return out;
}
inline bool ReadCsv(const std::string& path, Csv& c) {
    std::ifstream in(path.c_str());
    if (!in) { printf("smp::ReadCsv: cannot open %s\n", path.c_str()); return false; }
    std::string line; bool haveHead = false;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (!haveHead) { c.head = SplitCsv(line); haveHead = true; continue; }
        c.rows.push_back(SplitCsv(line));
    }
    return haveHead;
}

// ---- luminosity per day ----------------------------------------------------------
struct DayLumi {
    int day = 0, period = -1;
    double L6 = 0, L6stat = 0, L5 = 0;   // START-free, all runs of the day
    double LF2 = 0, LF3 = 0, LF32 = 0;   // PT2 / PT3 DAQ live fractions, ratio
    double S2fac6 = 1, f6 = 1;           // sector-2 deficit factor (inside L), day factor
    double Lstd_auth = 0, L5_auth = 0;   // published chain (LGAD START condition)
    // filled by ReadCoverage: runs of the day in the luminosity table
    int nRuns = 0, nRunsNoHF = 0;
};
inline std::map<int, DayLumi> ReadDayLumi() {
    std::map<int, DayLumi> m;
    Csv c;
    if (!ReadCsv(std::string(kLumiDir) + "/lumi_startfree_perday.csv", c)) return m;
    const int iD = c.col("day"), iL6 = c.col("L6"), iS = c.col("L6stat_data"), iL5 = c.col("L5");
    const int iLF2 = c.col("LF2"), iLF3 = c.col("LF3"), iR = c.col("LF3_over_LF2");
    const int iS2 = c.col("S2fac6"), iF6 = c.col("f6"), iA = c.col("Lstd_auth"), iA5 = c.col("L5_auth");
    for (const auto& r : c.rows) {
        DayLumi d;
        d.day = std::atoi(r[iD].c_str());
        d.period = PeriodOf(d.day);
        d.L6 = std::atof(r[iL6].c_str()); d.L6stat = std::atof(r[iS].c_str()); d.L5 = std::atof(r[iL5].c_str());
        d.LF2 = std::atof(r[iLF2].c_str()); d.LF3 = std::atof(r[iLF3].c_str()); d.LF32 = std::atof(r[iR].c_str());
        d.S2fac6 = std::atof(r[iS2].c_str()); d.f6 = std::atof(r[iF6].c_str());
        d.Lstd_auth = std::atof(r[iA].c_str()); d.L5_auth = std::atof(r[iA5].c_str());
        m[d.day] = d;
    }
    return m;
}

// ---- luminosity per period (sf_luminosity_perperiod.csv) ---------------------------
struct PerLumi {
    std::string name;
    double L6 = 0, L6stat = 0, L6unc = 0, L6hades = 0, L5 = 0, Lstd_auth = 0, L5_auth = 0;
};
inline std::vector<PerLumi> ReadPeriodLumi() {
    std::vector<PerLumi> v;
    Csv c;
    if (!ReadCsv(std::string(kLumiDir) + "/sf_luminosity_perperiod.csv", c)) return v;
    const int iP = c.col("period"), iL = c.col("L6_SF"), iS = c.col("L6_stat"), iU = c.col("L6_unc");
    const int iH = c.col("L6_hades"), i5 = c.col("L5_SF"), iA = c.col("Lstd_auth"), iA5 = c.col("L5_auth");
    for (const auto& r : c.rows) {
        PerLumi p;
        p.name = r[iP];
        p.L6 = std::atof(r[iL].c_str()); p.L6stat = std::atof(r[iS].c_str());
        p.L6unc = std::atof(r[iU].c_str()); p.L6hades = std::atof(r[iH].c_str());
        p.L5 = std::atof(r[i5].c_str()); p.Lstd_auth = std::atof(r[iA].c_str()); p.L5_auth = std::atof(r[iA5].c_str());
        v.push_back(p);
    }
    return v;
}

// ---- luminosity per run with the coverage of the FAT input -------------------------
struct RunLumi {
    std::string stem;
    int run = 0, day = 0, nList = 0, nMerged = 0, noHF = 0;
    double L6 = 0, L5 = 0, L6stat = 0;
};
// key = runnumber (= runid - 443000000, as stored in the PAT/FAT trees)
inline std::map<int, RunLumi> ReadRunLumi() {
    std::map<int, RunLumi> m;
    Csv c;
    if (!ReadCsv(std::string(kLumiDir) + "/lumi_coverage_perrun.csv", c)) return m;
    const int iS = c.col("dst_stem"), iR = c.col("runid"), iD = c.col("day"), iN = c.col("n_dst_list");
    const int iM = c.col("n_dst_merged"), iL6 = c.col("L6_SF"), iL5 = c.col("L5_SF"), iH = c.col("no_hf");
    for (const auto& r : c.rows) {
        RunLumi x;
        x.stem = r[iS];
        x.run = (int)(std::atol(r[iR].c_str()) - 443000000L);
        x.day = std::atoi(r[iD].c_str());
        x.nList = std::atoi(r[iN].c_str()); x.nMerged = std::atoi(r[iM].c_str());
        x.L6 = std::atof(r[iL6].c_str()); x.L5 = std::atof(r[iL5].c_str());
        x.noHF = std::atoi(r[iH].c_str());
        m[x.run] = x;
    }
    // statistical error per run from the per-run file
    Csv p;
    if (ReadCsv(std::string(kLumiDir) + "/lumi_startfree_perrun.csv", p)) {
        const int jR = p.col("runnumber"), jS = p.col("L6_stat");
        for (const auto& r : p.rows) {
            const int run = std::atoi(r[jR].c_str());
            auto it = m.find(run);
            if (it != m.end()) it->second.L6stat = std::atof(r[jS].c_str());
        }
    }
    return m;
}

// ---- runs present in the data (written by smp_pat_counts.C) --------------------------
// run -> rows of EpEm_ID, PipPimEpEm_ID (all, PT3), PipPimEpEp_ID (all, PT3), PipPimEmEm_ID (all, PT3)
struct RunRows { int day = 0; long long epem = 0, us0 = 0, us3 = 0, pp0 = 0, pp3 = 0, mm0 = 0, mm3 = 0; };
inline std::map<int, RunRows> ReadRunsPresent() {
    std::map<int, RunRows> m;
    for (int day : Days()) {
        std::ifstream in(Form("%s/pat/runs_%03d.txt", kWork, day));
        if (!in) { printf("smp::ReadRunsPresent: runs_%03d.txt missing (run smp_pat_counts.C)\n", day); continue; }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream is(line);
            int run; RunRows r; r.day = day;
            is >> run >> r.epem >> r.us0 >> r.us3 >> r.pp0 >> r.pp3 >> r.mm0 >> r.mm3;
            m[run] = r;
        }
    }
    return m;
}
// luminosity of the ANALYSED runs per day: sum of L6_SF over the runs that have
// EpEm_ID rows in the day file (17 557 runs; 6079.2 nb^-1). scaled = true scales
// a run by its fraction of merged DSTs (two runs: 6/10 and 8/10).
inline std::map<int, double> AnalysedLumiPerDay(bool scaled = false) {
    std::map<int, double> L;
    const std::map<int, RunLumi> rl = ReadRunLumi();
    const std::map<int, RunRows> rp = ReadRunsPresent();
    for (int day : Days()) L[day] = 0;
    for (const auto& p : rp) {
        if (p.second.epem <= 0) continue;
        auto it = rl.find(p.first);
        if (it == rl.end()) { printf("smp::AnalysedLumiPerDay: run %d not in the luminosity table\n", p.first); continue; }
        double w = 1.0;
        if (scaled && it->second.nList > 0) w = (double)it->second.nMerged / it->second.nList;
        L[it->second.day] += w * it->second.L6;
    }
    return L;
}

// ---- plot style and output ----------------------------------------------------------
inline void Style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42); gStyle->SetLabelFont(42, "xyz"); gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLegendFont(42);
    gStyle->SetLabelSize(0.050, "xyz"); gStyle->SetTitleSize(0.056, "xyz");
    gStyle->SetTitleOffset(1.05, "x"); gStyle->SetTitleOffset(1.25, "y");
    gStyle->SetPadLeftMargin(0.15); gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadTopMargin(0.08); gStyle->SetPadBottomMargin(0.14);
    gStyle->SetEndErrorSize(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);
}
// vector PDF cropped to the canvas (EPS -> ps2pdf -dEPSCrop) and a PNG preview; -dAutoRotatePages=/None:
// without it ghostscript turns a page whose text is mostly vertical by 90 deg (smp_cutflow; verification fix)
inline void Save(TCanvas* c, const std::string& base0) {
    const std::string n = base0.substr(base0.find_last_of('/') + 1);
    const bool inScope = std::string(kCntScope).find(" " + n + " ") != std::string::npos;
    const std::string eb = std::string(kCntScratch) + "/" + n;                         // EPS always in the scratch folder
    const std::string base = inScope ? std::string(kFigOut) + "/" + n + "_cnt" : eb;  // new name next to the original
    c->SaveAs((eb + ".eps").c_str());
    gSystem->Exec(Form("ps2pdf -dEPSCrop -dAutoRotatePages=/None %s.eps %s.pdf", eb.c_str(), base.c_str()));
    c->SaveAs((base + ".png").c_str());
}

} // namespace smp

#endif // SMP_COMMON_H
