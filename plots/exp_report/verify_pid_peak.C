// verify_pid_peak.C -- INDEPENDENT VERIFICATION of Table pid_peak_dstbeta (PID chapter): share of the eta signal whose
// leptons entered through the hadron pass of the DST track sorter.
//
// What it does
//   1. Loads the compact scan file with plots/ScanEngine.h (the standard spectra and CB machinery).
//   2. Takes the row class from verify_pid_class.root (written by verify_pid_rows.C from the PAT DST beta and the FAT
//      scan_nt rows, NOT from the chapter's tag trees): clsB = at least one lepton with DST beta < 0.9;
//      clsK = at least one lepton with kIsLepton == 0 (alternative definition with the sorter flag).
//      The class is stored in memory in Row::fold (fold = clsB + 2 clsK; fold is not used by the engine).
//   3. Selections written here from the definitions of the task: pre = vz > -160 (cor_ok == 1 is applied by the
//      engine; isBest >= 0, vz > -500, PT3 are FAT cuts); peak nominal = pre && oa_lab < 55 && M(pipi) > 0.30;
//      + OA(e+e-) > 2 deg; + vertex = both pions isUsedVertex; bump nominal = pre && oa_lab < 40 && 0.30 < M(pipi) < 0.41.
//   4. For all rows, class A and class B (each with its own like-sign CB = 2 sqrt(N++ N--), k = 1):
//      (a) P_main (ScanEngine FitPeakMain, frozen Gaussian 0.5439 / 5.4 MeV + pol2, 0.5225-0.6100, S = whole Gaussian);
//      (b) an independent simple estimate: pol2 fit of N - CB in the side bands 0.5225-0.5325 and 0.5550-0.6100
//          (TH1::Fit, weights 1/var), window sum of N - CB in 0.53364-0.55416 minus the integral of the pol2,
//          divided by the Gaussian content of the window (erf(1.9/sqrt2) = 0.9426);
//      (c) bump window W = N - CB in 0.4550-0.5225.
//      r = S_B / (S_A + S_B) with the error from the two independent fits; CB of both classes x 1.03 / 0.97.
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root,
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/verify/verify_pid_class.root
// Output : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/verify/verify_pid_peak.txt
// Run:
//   cd /home/przygoda/HADES/FAT/FAT/plots && root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_pid_verify", kTRUE)' \
//        'exp_report/verify_pid_peak.C+'
#include "../ScanEngine.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <cstdio>
#include <cmath>
#include <unordered_map>

namespace vpp {
const char* kCls = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/verify/verify_pid_class.root";
const char* kOut = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/verify/verify_pid_peak.txt";

// side-band estimate of the whole-Gaussian content (method b)
void SideBand(const scan::Spectra& s, double& S, double& eS)
{
   TH1D* h = s.Hist("sig", Form("hsb_%p", (const void*)&s));
   TH1D hs(*h);
   hs.SetName(Form("hsb2_%p", (const void*)&s));
   for (int i = 1; i <= hs.GetNbinsX(); ++i) {
      const double x = hs.GetBinCenter(i);
      const bool side = (x > 0.5225 && x < 0.5325) || (x > 0.5550 && x < 0.6100);
      if (!side) { hs.SetBinContent(i, 0); hs.SetBinError(i, 0); }
   }
   TF1 f("fsb", "pol2", 0.5225, 0.6100);
   hs.Fit(&f, "Q0 R");        // chi2 with the bin errors sqrt(N + var CB); bins outside the side bands are empty and skipped
   double win = 0, wvar = 0;
   for (int i = 0; i < s.nbin; ++i) {
      const double u = s.Frac(i, scan::kPeakLo, scan::kPeakHi);
      if (u <= 0) continue;
      win += u * s.sig[i];
      wvar += u * s.sigVar[i];
   }
   const double bg = f.Integral(scan::kPeakLo, scan::kPeakHi) / s.bw;
   const double gfrac = TMath::Erf(scan::kWinNsig / std::sqrt(2.));
   S = (win - bg) / gfrac;
   eS = std::sqrt(wvar) / gfrac;     // window statistics only (background-fit error neglected)
   delete h;
}

scan::Spectra ScaleCB(const scan::Spectra& in, double f)
{
   scan::Spectra s = in;
   for (int i = 0; i < s.nbin; ++i) { s.cb[i] *= f; s.cbVar[i] *= f * f; s.sig[i] = s.n[0][i] - s.cb[i]; s.sigVar[i] = s.n[0][i] + s.cbVar[i]; }
   for (auto& w : s.win) { w.cb *= f; w.cbVar *= f * f; }
   return s;
}
}

void verify_pid_peak()
{
   using namespace vpp;
   scan::Data d;
   if (!scan::Load(d, scan::kCompact)) return;
   FILE* fo = fopen(kOut, "w");
   auto out = [&](const char* s) { printf("%s", s); fprintf(fo, "%s", s); };

   // ---- classes
   TFile fc(kCls, "READ");
   const char* tn[3] = {"cls_us", "cls_pp", "cls_mm"};
   for (int k = 0; k < 3; ++k) {
      TTree* t = (TTree*)fc.Get(tn[k]);
      int fid, pe, b, kk;
      t->SetBranchAddress("file_id", &fid); t->SetBranchAddress("pat_entry", &pe);
      t->SetBranchAddress("clsB", &b); t->SetBranchAddress("clsK", &kk);
      std::unordered_map<long long, int> m;
      m.reserve(2 * t->GetEntries());
      for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); m[(long long)fid * 100000000LL + pe] = b + 2 * kk; }
      long long nm = 0, nB = 0, nK = 0, nmiss = 0;
      for (auto& r : d.s[k]) {
         auto it = m.find((long long)r.file_id * 100000000LL + r.pat_entry);
         if (it == m.end()) { r.fold = -1; ++nmiss; continue; }
         r.fold = it->second; ++nm; nB += (it->second & 1); nK += (it->second >> 1);
      }
      out(Form("%s: class rows %lld, compact rows %zu, matched %lld, missing %lld; clsB share %.4f, clsK share %.4f\n", tn[k],
               (long long)t->GetEntries(), d.s[k].size(), nm, nmiss, (double)nB / nm, (double)nK / nm));
   }

   using scan::Row;
   auto pre   = [](const Row& r) { return r.vz > -160.f; };
   auto peak  = [=](const Row& r) { return pre(r) && r.oa_lab < 55.f && r.mpipi > 0.30f; };
   auto peakO = [=](const Row& r) { return peak(r) && r.oa > 2.f; };
   auto peakV = [=](const Row& r) { return peak(r) && r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; };
   auto bump  = [=](const Row& r) { return pre(r) && r.oa_lab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f; };
   auto bumpV = [=](const Row& r) { return bump(r) && r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; };
   std::vector<scan::Sel> base = {pre, peak, peakO, peakV, bump, bumpV};
   const char* bn[6] = {"pre", "peak nominal", "peak nominal OA>2", "peak nominal + vertex", "bump nominal", "bump nominal + vertex"};
   // variants: 0 isBest>=0 / clsB, 1 isBest==1 / clsB, 2 isBest>=0 / clsK
   std::vector<scan::Sel> sels;
   for (int v = 0; v < 3; ++v)
      for (int b = 0; b < 6; ++b)
         for (int c = -1; c <= 1; ++c) {
            scan::Sel s = base[b];
            sels.push_back([s, c, v](const Row& r) {
               if (r.fold < 0) return false;
               if (v == 1 && r.isBest != 1) return false;
               const int cl = (v == 2) ? (r.fold >> 1) : (r.fold & 1);
               return (c < 0 || cl == c) && s(r);
            });
         }
   std::vector<scan::Spectra> sp = scan::FillMany(d, sels);
   const char* vn[3] = {"isBest >= 0, class by DST beta < 0.9", "isBest == 1, class by DST beta < 0.9", "isBest >= 0, class by kIsLepton == 0"};
   for (int v = 0; v < 3; ++v) {
      out(Form("\n=== %s\n", vn[v]));
      for (int b = 0; b < 6; ++b) {
         const bool isPeak = b < 4;
         double S[3], eS[3], Sb[3], eSb[3];
         for (int c = 0; c < 3; ++c) {
            const scan::Spectra& s = sp[(v * 6 + b) * 3 + c];
            if (isPeak) {
               const scan::PeakFit pf = scan::FitPeakMain(s);
               S[c] = pf.S; eS[c] = pf.eS;
               SideBand(s, Sb[c], eSb[c]);
               out(Form("  %-22s %-3s rows %8.0f  N(win) %6.0f CB(win) %7.1f  P_main S = %7.1f +- %5.1f  side-band S = %7.1f +- %5.1f\n", bn[b],
                        c == 0 ? "all" : (c == 1 ? "A" : "B"), s.nsel[0], s.win[0].n[0], s.win[0].cb, pf.S, pf.eS, Sb[c], eSb[c]));
            } else {
               S[c] = s.win[1].sig(); eS[c] = s.win[1].sigErr(); Sb[c] = S[c]; eSb[c] = eS[c];
               out(Form("  %-22s %-3s rows %8.0f  W = N - CB = %8.1f +- %5.1f (N %.0f, CB %.1f)\n", bn[b], c == 0 ? "all" : (c == 1 ? "A" : "B"),
                        s.nsel[0], S[c], eS[c], s.win[1].n[0], s.win[1].cb));
            }
         }
         const double r = S[2] / (S[1] + S[2]);
         const double er = std::hypot(S[1] * eS[2], S[2] * eS[1]) / std::pow(S[1] + S[2], 2);
         const double rb = Sb[2] / (Sb[1] + Sb[2]);
         const double erb = std::hypot(Sb[1] * eSb[2], Sb[2] * eSb[1]) / std::pow(Sb[1] + Sb[2], 2);
         double rv[2];
         for (int q = 0; q < 2; ++q) {
            const double f = q == 0 ? 1.03 : 0.97;
            double x[2];
            for (int c = 1; c <= 2; ++c) {
               const scan::Spectra s = ScaleCB(sp[(v * 6 + b) * 3 + c], f);
               x[c - 1] = isPeak ? scan::FitPeakMain(s).S : s.win[1].sig();
            }
            rv[q] = x[1] / (x[0] + x[1]);
         }
         out(Form("  %-22s r = S_B/(S_A+S_B) = %.3f +- %.3f (CB x1.03 %.3f, x0.97 %.3f); side-band method r = %.3f +- %.3f; A+B %.1f vs all %.1f\n",
                  bn[b], r, er, rv[0], rv[1], rb, erb, S[1] + S[2], S[0]));
      }
   }
   fclose(fo);
}
