// pid_peak_dstbeta.C -- PID chapter of the experimental-data report: how much of the eta signal has a lepton
// that entered the analysis through the HADRON pass of the DST track sorter (START-based DST beta < 0.9).
//
// What it does
//   1. Loads the compact scan file (trees us, pp, mm = e+e-, e+e+, e-e- hypotheses) with plots/ScanEngine.h.
//   2. Reads the tag trees tag_us / tag_pp / tag_mm of pid_final.root (written by pid_final_fill.C): for every FAT row
//      (file_id, pat_entry) the DST beta of its two lepton slots. Class of a row: A = both leptons with DST beta >= 0.9
//      (both could have entered through the lepton pass), B = at least one lepton with DST beta < 0.9 (that lepton
//      entered through the hadron pass, kIsLepton = 0). The class is kept in memory in the field Row::fold (not used here).
//   3. For the selections pre, peak nominal (with and without OA(e+e-) > 2 deg), peak nominal + vertex, bump nominal and
//      bump nominal + vertex, and for all rows / class A / class B: CB = 2 sqrt(N++ N--) (ScanEngine), eta peak P_main
//      (frozen Gaussian 0.5439 GeV/c2, sigma 5.4 MeV, pol2, 0.5225-0.6100; S = whole Gaussian) with the model spread
//      (48 variations, p > 0.05), and the bump window sum W = N - CB in 0.4550-0.5225 GeV/c2 (data window; W is not a yield).
//      Each class has its own like-sign CB (the like-sign rows are split by the same rule).
//   4. Share of class B in the signal: r = S_B / (S_A + S_B), stat. error from the independent fits of A and B;
//      the same with the CB of both classes scaled by 1.03 and 0.97 (CB sensitivity, k = 1).
//   Selections (cor_ok == 1 and isBest >= 0 always; FAT: vz > -500, PT3):
//     pre            vz > -160
//     peak nominal   pre && oa_lab < 55 && M(pipi) > 0.30
//     peak nom. OA>2 peak nominal && OA(e+e-) > 2 deg
//     peak + vertex  peak nominal && pip_isUsedVertex && pim_isUsedVertex
//     bump nominal   pre && oa_lab < 40 && 0.30 < M(pipi) < 0.41
//     bump nom. OA>2 bump nominal && OA(e+e-) > 2 deg (added by the verifier 2026-11-14)
//     bump + vertex  bump nominal && both pions in the vertex fit
//   and the same with isBest == 1.
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root (scan::kCompact),
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/pid_final.root
// Outputs: /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/pid_peak_dstbeta.{txt,tex}
// Run:
//   cd /home/przygoda/HADES/FAT/FAT/plots && root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_pid", kTRUE)' \
//        'exp_report/pid_peak_dstbeta.C+'
#include "../ScanEngine.h"
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
#include <cmath>
#include <unordered_map>
#include <vector>
#include <string>

namespace ppd {
const char* kFinal = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/pid_final.root";
const char* kTab   = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/";

inline long long key(int fid, int pe) { return (long long)fid * 100000000LL + pe; }

// scale the CB of a spectrum by f (and its variance by f^2) and recompute N - CB
scan::Spectra scaleCB(const scan::Spectra& in, double f)
{
   scan::Spectra s = in;
   for (int i = 0; i < s.nbin; ++i) {
      s.cb[i] *= f; s.cbVar[i] *= f * f;
      s.sig[i] = s.n[0][i] - s.cb[i]; s.sigVar[i] = s.n[0][i] + s.cbVar[i];
   }
   for (int j = 0; j < s.nw; ++j) { s.cbW[j] *= f; s.cbWVar[j] *= f * f; }
   for (auto& w : s.win) { w.cb *= f; w.cbVar *= f * f; }
   return s;
}
}

void pid_peak_dstbeta()
{
   using namespace ppd;
   scan::Data d;
   if (!scan::Load(d, scan::kCompact)) return;

   // ---- classes from the tag trees
   TFile ff(kFinal, "READ");
   const char* tn[3] = {"tag_us", "tag_pp", "tag_mm"};
   long long nTag[3] = {0, 0, 0}, nMatch[3] = {0, 0, 0}, nB[3] = {0, 0, 0};
   for (int k = 0; k < 3; ++k) {
      TTree* t = (TTree*)ff.Get(tn[k]);
      int fid, pe; float b1, b2;
      t->SetBranchAddress("file_id", &fid); t->SetBranchAddress("pat_entry", &pe);
      t->SetBranchAddress("b1", &b1); t->SetBranchAddress("b2", &b2);
      std::unordered_map<long long, int> cls;
      cls.reserve(t->GetEntries() * 2);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         cls[key(fid, pe)] = (b1 < 0.9f || b2 < 0.9f) ? 1 : 0;
      }
      nTag[k] = t->GetEntries();
      for (auto& r : d.s[k]) {
         auto it = cls.find(key(r.file_id, r.pat_entry));
         if (it == cls.end()) { r.fold = -1; continue; }
         r.fold = it->second;        // overwritten in memory: 0 = class A, 1 = class B
         ++nMatch[k];
         if (it->second == 1) ++nB[k];
      }
      printf("%s: tag rows %lld, compact rows %zu, matched %lld, class B %lld (%.4f)\n", tn[k], nTag[k], d.s[k].size(), nMatch[k], nB[k],
             nMatch[k] ? (double)nB[k] / nMatch[k] : 0.);
   }

   // ---- selections
   using scan::Row;
   auto pre   = [](const Row& r) { return r.vz > -160.f; };
   auto peak  = [pre](const Row& r) { return pre(r) && r.oa_lab < 55.f && r.mpipi > 0.30f; };
   auto peakO = [peak](const Row& r) { return peak(r) && r.oa > 2.f; };
   auto peakV = [peak](const Row& r) { return peak(r) && r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; };
   auto bump  = [pre](const Row& r) { return pre(r) && r.oa_lab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f; };
   auto bumpV = [bump](const Row& r) { return bump(r) && r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; };
   auto bumpO = [bump](const Row& r) { return bump(r) && r.oa > 2.f; };   // added by the verifier: OA > 2 deg variant also for the bump
   std::vector<scan::Sel> base = {pre, peak, peakO, peakV, bump, bumpO, bumpV};
   const int NB = 7;
   const char* bname[NB] = {"pre", "peak nominal", "peak nominal, OA>2", "peak nominal + vertex", "bump nominal", "bump nominal, OA>2", "bump nominal + vertex"};
   const char* btex[NB] = {"pre", "peak nominal", "peak nominal, OA$>2^\\circ$", "peak nominal + vertex", "bump nominal", "bump nominal, OA$>2^\\circ$", "bump nominal + vertex"};
   std::vector<scan::Sel> sels;
   for (int v = 0; v < 2; ++v)                 // v = 0: isBest >= 0, 1: isBest == 1
      for (int b = 0; b < NB; ++b)
         for (int c = -1; c <= 1; ++c) {      // -1 all, 0 class A, 1 class B
            scan::Sel s = base[b];
            sels.push_back([s, c, v](const Row& r) { return (v == 0 || r.isBest == 1) && (c < 0 || r.fold == c) && s(r); });
         }
   std::vector<scan::Spectra> sp = scan::FillMany(d, sels);

   FILE* f = fopen(Form("%spid_peak_dstbeta.txt", kTab), "w");
   FILE* x = fopen(Form("%spid_peak_dstbeta.tex", kTab), "w");
   fprintf(f, "# eta signal split by the DST (START) beta of the leptons: class A = both leptons DST beta >= 0.9, class B = at least one < 0.9\n");
   fprintf(f, "# (class B leptons entered through the hadron pass of the track sorter). P_main fit (frozen Gaussian, pol2, 0.5225-0.6100), S +- stat +- model spread;\n");
   fprintf(f, "# W = N - CB in the bump window 0.4550-0.5225 (not a yield). CB = 2 sqrt(N++ N--) per class. Rows: cor_ok == 1, FAT: vz > -500, PT3.\n");
   fprintf(f, "# class shares of the matched rows: us %.4f, pp %.4f, mm %.4f\n", (double)nB[0] / nMatch[0], (double)nB[1] / nMatch[1], (double)nB[2] / nMatch[2]);
   fprintf(x, "%% generated by plots/exp_report/pid_peak_dstbeta.C\n\\begin{table}[H]\n\\centering\n");
   fprintf(x, "\\caption{The $\\eta$ signal split by the START-based DST velocity of the two leptons: class A, both leptons with "
              "$\\beta_{\\mathrm{DST}}\\geq 0.9$; class B, at least one lepton with $\\beta_{\\mathrm{DST}}<0.9$ (it entered through the hadron pass of the "
              "track sorter). Peak: P\\_main fit of $N-\\mathrm{CB}$ (frozen Gaussian, pol2), $S\\pm$stat$\\pm$model spread (n/a: fewer than two of the 48 "
              "model variations have $p>0.05$); share of class B "
              "$r=S_B/(S_A+S_B)$ with its statistical error and, in brackets, $r$ with the CB of both classes scaled by $\\pm3\\,\\%%$. "
              "Bump: $W=N-\\mathrm{CB}$ in the data window 0.4550--0.5225\\,GeV/$c^2$ (a visibility measure, not a yield). "
              "CB $=2\\sqrt{N_{++}N_{--}}$ separately for each class; isBest $\\geq 0$; cor\\_ok $=1$, $z>-160$\\,mm, PT3.}\n");
   fprintf(x, "\\label{tab:pid_peak_dstbeta}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{@{}l r r r l@{}}\n\\toprule\n");
   fprintf(x, "Selection & all & class A & class B & $r$ (class B share) \\\\\n\\midrule\n");
   for (int v = 0; v < 2; ++v) {
      fprintf(f, "\n=== %s\n", v == 0 ? "isBest >= 0" : "isBest == 1");
      for (int b = 0; b < NB; ++b) {
         const bool isPeak = (b < 4);
         double S[3], eS[3], sd[3];
         int np[3] = {0, 0, 0};
         for (int c = 0; c < 3; ++c) {
            const scan::Spectra& s = sp[(v * NB + b) * 3 + c];
            if (isPeak) {
               const scan::PeakFit pf = scan::FitPeakMain(s);
               const scan::Spread spr = scan::PeakModelSpread(s);
               S[c] = pf.S; eS[c] = pf.eS; sd[c] = spr.half; np[c] = spr.npass;
               fprintf(f, "%-24s %-5s S = %7.1f +- %5.1f +- %5.1f (spread from %2d of 48 fits with p > 0.05)  chi2/ndf %.1f/%d  N(win) %.0f CB(win) %.1f (CB/N %.3f)  rows %.0f / %.0f / %.0f\n", bname[b],
                       c == 0 ? "all" : (c == 1 ? "A" : "B"), pf.S, pf.eS, spr.half, spr.npass, pf.chi2, pf.ndf, pf.N, pf.CB, pf.N > 0 ? pf.CB / pf.N : 0., s.nsel[0], s.nsel[1], s.nsel[2]);
            } else {
               const scan::WinSum& w = s.win[1];
               S[c] = w.sig(); eS[c] = w.sigErr(); sd[c] = 0;
               fprintf(f, "%-24s %-5s W = %8.1f +- %6.1f  N %.0f CB %.1f (CB/N %.3f)  rows %.0f / %.0f / %.0f\n", bname[b], c == 0 ? "all" : (c == 1 ? "A" : "B"),
                       w.sig(), w.sigErr(), w.n[0], w.cb, w.cbShare(), s.nsel[0], s.nsel[1], s.nsel[2]);
            }
         }
         const double sum = S[1] + S[2];
         const double r = S[2] / sum;
         const double er = std::sqrt(std::pow(S[1] * eS[2], 2) + std::pow(S[2] * eS[1], 2)) / (sum * sum);
         // CB +- 3 %
         double rv[2];
         for (int q = 0; q < 2; ++q) {
            const double fc = q == 0 ? 1.03 : 0.97;
            double SS[2];
            for (int c = 1; c <= 2; ++c) {
               const scan::Spectra s = scaleCB(sp[(v * NB + b) * 3 + c], fc);
               SS[c - 1] = isPeak ? scan::FitPeakMain(s).S : s.win[1].sig();
            }
            rv[q] = SS[1] / (SS[0] + SS[1]);
         }
         fprintf(f, "%-24s share of class B: r = %.3f +- %.3f (CB x1.03: %.3f, x0.97: %.3f); A + B = %.1f vs all %.1f\n", bname[b], r, er, rv[0], rv[1], sum, S[0]);
         if (v == 0) {
            if (isPeak)
            {
               TString c3[3];
               for (int c = 0; c < 3; ++c)
                  c3[c] = np[c] >= 2 ? Form("%.0f $\\pm$ %.0f $\\pm$ %.0f", S[c], eS[c], sd[c]) : Form("%.0f $\\pm$ %.0f $\\pm$ n/a", S[c], eS[c]);
               fprintf(x, "%s & %s & %s & %s & %.3f $\\pm$ %.3f (%.3f, %.3f) \\\\\n", btex[b], c3[0].Data(), c3[1].Data(), c3[2].Data(), r, er, rv[0], rv[1]);
            }
            else
               fprintf(x, "%s ($W$) & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.3f $\\pm$ %.3f (%.3f, %.3f) \\\\\n", btex[b], S[0], eS[0], S[1], eS[1],
                       S[2], eS[2], r, er, rv[0], rv[1]);
         }
      }
   }
   fprintf(x, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
   fclose(f);
   fclose(x);
   printf("written %spid_peak_dstbeta.{txt,tex}\n", kTab);
}
