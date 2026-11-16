// spec_vz500_cnt.C -- COPY of plots/exp_report/spec_vz500.C made by the final checker (integ4, 2026-11-15): figure strings
// 'rows' -> 'counts' (analyst's terminology); figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png};
// tables, ROOT files and EPS -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/check/figs/spec. Physics code unchanged.
// ============================================================================
// spec_vz500.C  --  effect of the mandatory FAT cut "event vertex z > -500 mm"
//                   (and of PT3, cor_ok, vz > -160 mm) on unlike-sign versus
//                   like-sign rows, measured on the PAT input (before FAT)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// The FAT production (and therefore the compact file of all other spec_*
// macros) contains only rows with eVertReco_z > -500 mm and PT3; this macro
// reads the 29 PAT day files directly to show what these cuts remove.
//
// What it does
//   mode "fill" (part k of n): for the files with index i % n == k of
//     /hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list it reads the trees
//     PipPimEpEm_ID (us), PipPimEpEp_ID (pp), PipPimEmEm_ID (mm), only the
//     branches needed (isBest, trigbit, runnumber, event, eVertReco_z, the
//     four momenta p_corr / theta / phi as FAT reads them, the pion
//     isUsedVertex flags), computes M4 = M(pi+pi-e+e-) with the energy-loss
//     corrected momenta exactly as FAT (lepton slot 1: <p1>_p_corr_ep, slot 2:
//     <p2>_p_corr_em, pions <pi>_p_corr_<pi>; cor_ok = all four finite) and
//     counts rows at the stages
//        0: all PAT rows (isBest >= 0)
//        1: + eVertReco_z > -500 mm               (FAT)
//        2: + PT3, trigbit == 8192                (FAT)
//        3: + cor_ok                              (analysis)
//        4: + vz > -160 mm                        (analysis: = pre)
//     in total, in the standard M4 regions (cor_ok rows) and in 10 MeV M4 bins,
//     with event groups = adjacent rows with equal (run, event) inside one file
//     (ratio-estimator errors of the pass fractions); also the vz distribution
//     of all rows (-1000 ... +200 mm). Writes report_exp/spec/vz500/part_<k>.root.
//   mode "draw": adds the parts, writes the table and the figure.
//
// Inputs : /hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list (29 PAT day files, 374 GB)
// Outputs: /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/spec/vz500/part_<k>.root
//          /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/spec_vz500.txt, spec_cutflow.tex
//          /home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/spec_vz500.pdf (+ .png)
//
// Run (from the repository root /home/przygoda/HADES/FAT/FAT): six fill jobs in parallel, then draw
//   for k in 0 1 2 3 4 5; do root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_spec", kTRUE)' "plots/exp_report/spec_vz500.C+(\"fill\",$k,6)" > /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/spec/vz500/fill_$k.log 2>&1 & done; wait
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_spec", kTRUE)' 'plots/exp_report/spec_vz500.C+("draw")'
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

#include "spec_common_cnt2.h"

namespace specvz {

using namespace spec;

const char* const kList = "/hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list";
const std::string kDir = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/spec/vz500";   // original part files, read only (draw mode)
const int kNStage = 5;
const char* const kStage[kNStage] = {"PAT rows (isBest >= 0)", "+ vertex z > -500 mm (FAT)", "+ PT3 trigbit == 8192 (FAT)",
                                     "+ cor_ok (analysis)", "+ vz > -160 mm (= pre)"};
const char* const kTree[3] = {"PipPimEpEm_ID", "PipPimEpEp_ID", "PipPimEmEm_ID"};
const char* const kPre1[3] = {"ep", "ep1", "em1"};
const char* const kPre2[3] = {"em", "ep2", "em2"};
// M4 regions as in spec_m4.C (cor_ok rows)
const int kNReg = 9;
const double kRegLo[kNReg] = {0.0, 0.455, 0.5225, scan::kPeakLo, scan::kPeakHi, 0.56, 0.70, 0.90, 1.10};
const double kRegHi[kNReg] = {0.455, 0.5225, scan::kPeakLo, scan::kPeakHi, 0.56, 0.70, 0.90, 1.10, 10.0};
const char* const kRegName[kNReg] = {"M4 < 0.455", "bump window", "gap", "peak window", "upper gap", "0.56-0.70",
                                     "0.70-0.90", "0.90-1.10", "M4 > 1.10"};
const double kMlo = 0.35, kMhi = 1.10; const int kMnb = 75;    // 10 MeV bins for the figure
const double kMe = 0.51099895, kMpi = 139.57039;               // MeV/c2 (FAT pparticle.h)

// accumulators of one sample: per bin (region or M4 bin) and transition (stage s -> s+1):
// sum kp, kt, kp^2, kt^2, kp*kt over event groups (ratio-estimator variance)
struct Acc {
    TH1D* h[5];   // p, t, pp, tt, pt
    void Make(const char* name, int nb, double lo, double hi) {
        const char* c[5] = {"p", "t", "pp", "tt", "pt"};
        for (int i = 0; i < 5; ++i) { h[i] = new TH1D(Form("%s_%s", name, c[i]), "", nb, lo, hi); h[i]->SetDirectory(nullptr); }
    }
    void Load(TFile* f, const char* name, bool add) {
        const char* c[5] = {"p", "t", "pp", "tt", "pt"};
        for (int i = 0; i < 5; ++i) {
            TH1D* x = (TH1D*)f->Get(Form("%s_%s", name, c[i]));
            if (add) h[i]->Add(x);
            else { h[i] = (TH1D*)x->Clone(); h[i]->SetDirectory(nullptr); }
        }
    }
    void Write() { for (int i = 0; i < 5; ++i) h[i]->Write(); }
    double F(int b) const { const double t = h[1]->GetBinContent(b + 1); return t > 0 ? h[0]->GetBinContent(b + 1) / t : 0; }
    double E(int b) const {
        const double t = h[1]->GetBinContent(b + 1);
        if (t <= 0) return 0;
        const double f = F(b);
        const double v = h[2]->GetBinContent(b + 1) - 2 * f * h[4]->GetBinContent(b + 1) + f * f * h[3]->GetBinContent(b + 1);
        return std::sqrt(std::max(v, 0.)) / t;
    }
    double T(int b) const { return h[1]->GetBinContent(b + 1); }
    double P(int b) const { return h[0]->GetBinContent(b + 1); }
};

std::vector<std::string> ReadList() {
    std::vector<std::string> v; std::ifstream in(kList); std::string s;
    while (std::getline(in, s)) if (!s.empty() && s[0] != '#') v.push_back(s);
    return v;
}

void Fill(int part, int nparts) {
    gSystem->mkdir(kDir.c_str(), kTRUE);
    const std::vector<std::string> files = ReadList();
    TH1D* hCount[3]; TH1D* hVz[3];
    Acc aReg[3][kNStage - 1], aM4[3][kNStage - 1];
    for (int k = 0; k < 3; ++k) {
        hCount[k] = new TH1D(Form("count_%d", k), "", kNStage, 0, kNStage); hCount[k]->SetDirectory(nullptr);
        hVz[k] = new TH1D(Form("vz_%d", k), "", 240, -1000, 200); hVz[k]->SetDirectory(nullptr);
        for (int s = 0; s + 1 < kNStage; ++s) {
            aReg[k][s].Make(Form("reg_%d_%d", k, s), kNReg, 0, kNReg);
            aM4[k][s].Make(Form("m4_%d_%d", k, s), kMnb, kMlo, kMhi);
        }
    }
    TStopwatch sw;
    for (size_t fi = 0; fi < files.size(); ++fi) {
        if ((int)(fi % nparts) != part) continue;
        TFile* f = TFile::Open(files[fi].c_str(), "READ");
        if (!f || f->IsZombie()) { printf("cannot open %s\n", files[fi].c_str()); continue; }
        for (int k = 0; k < 3; ++k) {
            TTree* t = (TTree*)f->Get(kTree[k]);
            if (!t) { printf("tree %s missing in %s\n", kTree[k], files[fi].c_str()); continue; }
            t->SetBranchStatus("*", 0);
            std::vector<std::string> names = {"isBest", "trigbit", "runnumber", "event", "eVertReco_z",
                "pip_p_corr_pip", "pip_theta", "pip_phi", "pim_p_corr_pim", "pim_theta", "pim_phi",
                std::string(kPre1[k]) + "_p_corr_ep", std::string(kPre1[k]) + "_theta", std::string(kPre1[k]) + "_phi",
                std::string(kPre2[k]) + "_p_corr_em", std::string(kPre2[k]) + "_theta", std::string(kPre2[k]) + "_phi"};
            std::vector<Float_t> val(names.size(), 0.f);
            bool ok = true;
            for (size_t i = 0; i < names.size(); ++i) {
                if (!t->GetBranch(names[i].c_str())) { printf("branch %s missing in %s\n", names[i].c_str(), kTree[k]); ok = false; continue; }
                t->SetBranchStatus(names[i].c_str(), 1);
                t->SetBranchAddress(names[i].c_str(), &val[i]);
            }
            if (!ok) continue;
            const Long64_t n = t->GetEntries();
            // group bookkeeping: (run, event) of the previous row; per group a list of (stage transition, bin kind, bin, pass)
            double curRun = -1, curEvt = -1;
            std::vector<std::array<int, 4>> g;   // transition s, kind (0 region, 1 M4 bin), bin, pass
            auto flush = [&]() {
                // aggregate kp, kt per (s, kind, bin)
                std::sort(g.begin(), g.end());
                size_t i = 0;
                while (i < g.size()) {
                    size_t j = i; int kp = 0, kt = 0;
                    while (j < g.size() && g[j][0] == g[i][0] && g[j][1] == g[i][1] && g[j][2] == g[i][2]) { kp += g[j][3]; kt += 1; ++j; }
                    Acc& a = (g[i][1] == 0) ? aReg[k][g[i][0]] : aM4[k][g[i][0]];
                    const int b = g[i][2] + 1;
                    a.h[0]->AddBinContent(b, kp); a.h[1]->AddBinContent(b, kt);
                    a.h[2]->AddBinContent(b, (double)kp * kp); a.h[3]->AddBinContent(b, (double)kt * kt); a.h[4]->AddBinContent(b, (double)kp * kt);
                    i = j;
                }
                g.clear();
            };
            for (Long64_t e = 0; e < n; ++e) {
                t->GetEntry(e);
                const double run = val[2], evt = val[3];
                if (run != curRun || evt != curEvt) { flush(); curRun = run; curEvt = evt; }
                if (!(val[0] > -0.5)) continue;                        // isBest >= 0
                const double vz = val[4];
                hVz[k]->Fill(std::max(-999.9, std::min(199.9, vz)));
                const bool st[kNStage] = {true, vz > -500., vz > -500. && std::lround(val[1]) == 8192, false, false};
                // corrected four-momenta as in FAT (MeV/c, deg)
                double px = 0, py = 0, pz = 0, E = 0; bool cor = true;
                const int ip[4] = {5, 8, 11, 14};
                const double ms[4] = {kMpi, kMpi, kMe, kMe};
                for (int q = 0; q < 4; ++q) {
                    const double p = val[ip[q]], th = val[ip[q] + 1] * TMath::DegToRad(), ph = val[ip[q] + 2] * TMath::DegToRad();
                    if (!std::isfinite(p)) { cor = false; break; }
                    px += p * std::sin(th) * std::cos(ph); py += p * std::sin(th) * std::sin(ph); pz += p * std::cos(th);
                    E += std::sqrt(p * p + ms[q] * ms[q]);
                }
                const double m4 = cor ? 1e-3 * std::sqrt(std::max(E * E - px * px - py * py - pz * pz, 0.)) : -1;
                bool s3 = st[2] && cor, s4 = s3 && vz > -160.;
                const bool stage[kNStage] = {st[0], st[1], st[2], s3, s4};
                for (int s = 0; s < kNStage; ++s) if (stage[s]) hCount[k]->Fill(s + 0.5);
                // pass fractions of the transitions s -> s+1 (only defined with an M4: cor_ok rows; the
                // transition 2 -> 3 is cor_ok itself and is counted in the region of the REC-free row: skipped)
                if (!cor) continue;
                int reg = -1;
                for (int q = 0; q < kNReg; ++q) if (m4 >= kRegLo[q] && m4 < kRegHi[q]) { reg = q; break; }
                const int mb = (m4 >= kMlo && m4 < kMhi) ? (int)((m4 - kMlo) / (kMhi - kMlo) * kMnb) : -1;
                const bool cs[kNStage] = {true, st[1], st[2], st[2], s4};   // for cor_ok rows stage 3 == stage 2
                for (int s = 0; s + 1 < kNStage; ++s) {
                    if (!cs[s]) continue;
                    const int pass = cs[s + 1] ? 1 : 0;
                    if (reg >= 0) g.push_back({s, 0, reg, pass});
                    if (mb >= 0) g.push_back({s, 1, mb, pass});
                }
            }
            flush();
            printf("file %zu %s tree %s: %lld rows, %.0f s\n", fi, files[fi].c_str(), kTree[k], n, sw.RealTime());
            sw.Continue();
            t->ResetBranchAddresses();
        }
        f->Close(); delete f;
    }
    TFile* o = TFile::Open(Form("%s/part_%d.root", kDir.c_str(), part), "RECREATE");
    for (int k = 0; k < 3; ++k) {
        hCount[k]->Write(); hVz[k]->Write();
        for (int s = 0; s + 1 < kNStage; ++s) { aReg[k][s].Write(); aM4[k][s].Write(); }
    }
    o->Close();
    printf("part %d of %d done in %.0f s\n", part, nparts, sw.RealTime());
}

FILE* gOut = nullptr;
void Pr(const char* fmt, ...) {
    va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
    if (gOut) { va_start(a, fmt); vfprintf(gOut, fmt, a); va_end(a); }
}

void Draw() {
    Style();
    TH1D* hCount[3] = {nullptr, nullptr, nullptr}; TH1D* hVz[3] = {nullptr, nullptr, nullptr};
    Acc aReg[3][kNStage - 1], aM4[3][kNStage - 1];
    int nparts = 0;
    for (int part = 0; part < 64; ++part) {
        const TString fn = Form("%s/part_%d.root", kDir.c_str(), part);
        if (gSystem->AccessPathName(fn)) continue;
        TFile* f = TFile::Open(fn);
        const bool add = (nparts > 0);
        for (int k = 0; k < 3; ++k) {
            TH1D* c = (TH1D*)f->Get(Form("count_%d", k)); TH1D* v = (TH1D*)f->Get(Form("vz_%d", k));
            if (add) { hCount[k]->Add(c); hVz[k]->Add(v); }
            else { hCount[k] = (TH1D*)c->Clone(); hCount[k]->SetDirectory(nullptr); hVz[k] = (TH1D*)v->Clone(); hVz[k]->SetDirectory(nullptr); }
            for (int s = 0; s + 1 < kNStage; ++s) { aReg[k][s].Load(f, Form("reg_%d_%d", k, s), add); aM4[k][s].Load(f, Form("m4_%d_%d", k, s), add); }
        }
        f->Close(); delete f;
        nparts++;
    }
    if (!nparts) { printf("no parts in %s\n", kDir.c_str()); return; }

    gOut = fopen(Form("%s/spec_vz500.txt", kTabOut), "w");
    Pr("# spec_vz500.txt -- what the FAT preselection (vertex z > -500 mm, PT3) and the analysis cuts (cor_ok, vz > -160 mm) remove,\n");
    Pr("#   unlike-sign (us, PipPimEpEm_ID) vs like-sign (pp, PipPimEpEp_ID; mm, PipPimEmEm_ID) rows; written by plots/exp_report/spec_vz500.C\n");
    Pr("# input: the 29 PAT day files of %s (%d parts added); all rows have isBest 0 or 1\n", kList, nparts);
    Pr("# A. cut flow (rows)\n");
    Pr("   %-30s %12s %8s %10s %8s %10s %8s %8s\n", "stage", "us", "kept", "pp", "kept", "mm", "kept", "CB/N*");
    for (int s = 0; s < kNStage; ++s) {
        const double u = hCount[0]->GetBinContent(s + 1), p = hCount[1]->GetBinContent(s + 1), m = hCount[2]->GetBinContent(s + 1);
        const double u0 = s ? hCount[0]->GetBinContent(s) : u, p0 = s ? hCount[1]->GetBinContent(s) : p, m0 = s ? hCount[2]->GetBinContent(s) : m;
        Pr("   %-30s %12.0f %8.4f %10.0f %8.4f %10.0f %8.4f %8.4f\n", kStage[s], u, u / u0, p, p / p0, m, m / m0, 2 * std::sqrt(p * m) / u);
    }
    Pr("   * CB/N = 2 sqrt(N++ N--)/N of the totals (not the M4-binned CB); 'kept' = fraction kept by this stage.\n");
    Pr("\n# B. row pass fractions per M4 region (corrected momenta, cor_ok rows), errors with the event-group correlation\n");
    const char* tr[kNStage - 1] = {"vertex z > -500 mm (base: PAT rows)", "PT3 trigbit == 8192 (base: z > -500)",
                                   "cor_ok (== 1 here by construction)", "vz > -160 mm (base: z > -500, PT3, cor_ok)"};
    for (int s = 0; s + 1 < kNStage; ++s) {
        if (s == 2) continue;
        Pr("== %s\n", tr[s]);
        Pr("   %-14s %20s %20s %20s %9s %8s\n", "region", "f(us)", "f(pp)", "f(mm)", "CBkept", "CB/N chg");
        for (int r = 0; r < kNReg; ++r) {
            const double cbk = std::sqrt(aReg[1][s].F(r) * aReg[2][s].F(r));
            Pr("   %-14s %12.4f +- %5.4f %12.4f +- %5.4f %12.4f +- %5.4f %9.4f %8.3f\n", kRegName[r], aReg[0][s].F(r), aReg[0][s].E(r),
               aReg[1][s].F(r), aReg[1][s].E(r), aReg[2][s].F(r), aReg[2][s].E(r), cbk, aReg[0][s].F(r) > 0 ? cbk / aReg[0][s].F(r) : 0.);
        }
    }
    Pr("   CBkept = sqrt(f(pp) f(mm)) = kept fraction of 2 sqrt(N++ N--); 'CB/N chg' = factor by which the cut changes CB/N.\n");
    // vz below -500
    for (int k = 0; k < 3; ++k) {
        const int b500 = hVz[k]->FindBin(-500.001);
        Pr("# C. vz of all PAT rows, sample %s: fraction at vz <= -500 mm %.4f (of which at the lower clamp -1000: %.4f)\n",
           k == 0 ? "us" : (k == 1 ? "pp" : "mm"), hVz[k]->Integral(1, b500) / hVz[k]->Integral(), hVz[k]->GetBinContent(1) / hVz[k]->Integral());
    }
    fclose(gOut); gOut = nullptr;
    // LaTeX cut flow
    {
        FILE* o = fopen(Form("%s/spec_cutflow.tex", kTabOut), "w");
        fprintf(o, "%% written by plots/exp_report/spec_vz500.C -- row cut flow from the PAT input to pre, three hypotheses\n");
        fprintf(o, "\\begin{table}[H]\n\\centering\n\\scriptsize\n");
        fprintf(o, "\\caption{Rows (stored particle combinations, $\\mathrm{isBest}\\ge0$) of the three PAT hypotheses of the 29 day files, "
                   "from the PAT input to the selection \\emph{pre}. In brackets: the fraction kept by the stage. The event-vertex cut "
                   "$z>-500$~mm removes the rows without a fitted event vertex (PAT default $z=-1000$~mm). CB/$N$ $=2\\sqrt{N_{++}N_{--}}/N$ "
                   "of the totals (not M4-binned).}\n");
        fprintf(o, "\\label{tab:spec_cutflow}\n\\begin{tabular}{lrrrr}\n\\toprule\n");
        fprintf(o, "stage & $e^+e^-$ & $e^+e^+$ & $e^-e^-$ & CB/$N$ \\\\\n\\midrule\n");
        const char* st[kNStage] = {"PAT rows", "+ event vertex $z>-500$~mm (FAT)", "+ PT3, trigbit $=8192$ (FAT)",
                                   "+ $\\mathrm{cor\\_ok}=1$", "+ $v_z>-160$~mm (= \\emph{pre})"};
        for (int s = 0; s < kNStage; ++s) {
            double v[3], f[3];
            for (int k = 0; k < 3; ++k) { v[k] = hCount[k]->GetBinContent(s + 1); f[k] = s ? v[k] / hCount[k]->GetBinContent(s) : 1.0; }
            fprintf(o, "%s & %.0f (%.3f) & %.0f (%.3f) & %.0f (%.3f) & %.4f \\\\\n", st[s], v[0], f[0], v[1], f[1], v[2], f[2], 2 * std::sqrt(v[1] * v[2]) / v[0]);
        }
        fprintf(o, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
        fclose(o);
    }

    // ---- figure: (a) vz of all PAT rows, (b) pass fraction of z > -500 versus M4 ----
    TCanvas* c = new TCanvas("c_vz500", "", 1800, 700);
    c->Divide(2, 1, 0.001, 0.001);
    const int col[3] = {kBlack, kGreen + 2, kOrange + 7};   // [review fix] like-sign key of the report: e+e+ green squares, e-e- orange triangles
    const int mk[3] = {20, 25, 26};
    const char* smp[3] = {"e^{+}e^{-} (unlike-sign)", "e^{+}e^{+}", "e^{-}e^{-}"};
    c->cd(1); gPad->SetLogy();
    double ym = 0, fNoVtx[3];
    for (int k = 0; k < 3; ++k) {
        // first bin = rows with vz <= -1000 (no fitted event vertex: eVertReco_z = -1000 default); quoted as a number
        fNoVtx[k] = hVz[k]->GetBinContent(1) / hVz[k]->Integral();
        hVz[k]->Scale(1.0 / hVz[k]->Integral(), "width");     // fraction of ALL rows of the sample, per mm
        hVz[k]->SetBinContent(1, 0); hVz[k]->SetBinError(1, 0);
        // last bin = all rows with vz >= 195 mm (the fill clamps them to 199.9): not a density, not drawn
        // (verifier, 2026-11-14)
        hVz[k]->SetBinContent(hVz[k]->GetNbinsX(), 0); hVz[k]->SetBinError(hVz[k]->GetNbinsX(), 0);
        Deco(hVz[k], col[k], mk[k], 0.6);
        ym = std::max(ym, hVz[k]->GetMaximum());
    }
    TH1D* fr = new TH1D("fr_a", "", 100, -995, 200); fr->SetDirectory(nullptr);
    fr->SetMinimum(1e-7); fr->SetMaximum(5000 * ym);
    fr->GetXaxis()->SetTitle(Fx("event vertex v_{z} [mm]").c_str()); fr->GetYaxis()->SetTitle("fraction of counts / mm");   // Scale(..., "width") gives per mm (verifier fix; was "/ 5 mm")
    fr->GetXaxis()->SetNdivisions(506);
    fr->Draw("axis");
    for (int k = 0; k < 3; ++k) hVz[k]->Draw("e same");
    VLine(-500, 1e-7, 1.2 * ym, kRed + 1, 2); VLine(-160, 1e-7, 1.2 * ym, kGray + 2, 7);
    Txt(0.16, 0.935, "(a) PAT counts, before the FAT preselection", 0.050);
    TLegend* lg = Leg(0.18, 0.74, 0.58, 0.90, 0.042);
    for (int k = 0; k < 3; ++k) LE(lg, hVz[k], smp[k], "pe");
    lg->Draw();
    Txt(0.18, 0.685, "not shown: no fitted event vertex (v_{z} = #minus1000 mm)", 0.040);
    Txt(0.18, 0.635, Form("e^{+}e^{-} %.1f %%, e^{+}e^{+} %.1f %%, e^{-}e^{-} %.1f %% of the counts", 100 * fNoVtx[0], 100 * fNoVtx[1], 100 * fNoVtx[2]), 0.040);
    Txt(0.485, 0.50, "keep v_{z} > #minus500 #rightarrow", 0.038, 31, kRed + 1);
    c->cd(2);
    TH1D* fr2 = new TH1D("fr_b", "", 100, kMlo, kMhi); fr2->SetDirectory(nullptr);
    fr2->SetMinimum(0); fr2->SetMaximum(1.25);
    fr2->GetXaxis()->SetTitle(Fx("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]").c_str());
    fr2->GetYaxis()->SetTitle("pass fraction of the counts");
    fr2->GetXaxis()->SetNdivisions(506);
    fr2->Draw("axis");
    Band(scan::kBumpLo, scan::kBumpHi, 0, 1, ColBump()); Band(scan::kPeakLo, scan::kPeakHi, 0, 1, ColPeak());
    TGraphErrors* gk[3];
    for (int k = 0; k < 3; ++k) {
        TGraphErrors* g = new TGraphErrors();
        g->SetMarkerColor(col[k]); g->SetLineColor(col[k]); g->SetMarkerStyle(mk[k]); g->SetMarkerSize(0.8);
        for (int b = 0; b < kMnb; ++b) {
            if (aM4[k][0].T(b) < 5) continue;
            g->SetPoint(g->GetN(), kMlo + (b + 0.5) * (kMhi - kMlo) / kMnb, aM4[k][0].F(b));
            g->SetPointError(g->GetN() - 1, 0, aM4[k][0].E(b));
        }
        g->Draw("p same"); gk[k] = g;
    }
    gPad->RedrawAxis();
    Txt(0.16, 0.935, "(b) v_{z} > #minus500 mm versus M_{4} (corrected momenta)", 0.050);
    TLegend* l2 = Leg(0.18, 0.80, 0.96, 0.89, 0.042); l2->SetNColumns(3);
    for (int k = 0; k < 3; ++k) LE(l2, gk[k], smp[k], "pe");
    l2->Draw();
    Save(c, "spec_vz500");
}

} // namespace specvz

void spec_vz500_cnt(const char* mode = "draw", int part = 0, int nparts = 1) {
    const std::string m = mode;
    if (m == "fill") { printf("spec_vz500_cnt: draw mode only\n"); return; }
    else specvz::Draw();
}
