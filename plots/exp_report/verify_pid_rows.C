// verify_pid_rows.C -- INDEPENDENT VERIFICATION of the PID chapter (prefix pid) of the experimental-data report.
//
// What it does (written independently of pid_final_fill.C / pid_draw.C; different input path for the selection):
//   1. Reads the START-based DST beta (HParticleCand::getBeta, PAT column <slot>_beta) and the DST tof of the lepton
//      and pion slots directly from the 29 PAT day files (trees PipPimEpEm_ID, PipPimEpEp_ID, PipPimEmEm_ID).
//   2. Takes the row selection from the FAT outputs (tree scan_nt of output_pippim{epem,epep,emem}_exp.root; NOT
//      from the compact scan file): cor_ok, eVertReco_z, isBest, isUsedVertex of the pions, the corrected masses
//      m_pippimepem_cor, m_pippim_cor and oa_pippim_epem_lab_cor, and p_rec / beta_new of the four slots; the PAT
//      row is found by (file_id = day, pat_entry).
//   3. Selections: fat = all scan_nt rows (isBest >= 0, vz > -500, PT3 are FAT cuts, checked: 0 violations);
//      pre = fat && cor_ok == 1 && vz > -160; pre_best = pre && isBest == 1; pre_vtx = pre && both pions
//      isUsedVertex; peak = pre && oa_lab < 55 && M(pipi) > 0.30 && 0.53364 <= M4 < 0.55416 (data peak window);
//      bump = pre && oa_lab < 40 && 0.30 < M(pipi) < 0.41 && 0.4550 <= M4 < 0.5225 (data bump window);
//      peak_oa2, bump_oa2 = the same with OA(e+e-) > 2 deg.
//   4. Per selection: rows, share of e+ / e- with DST beta < 0.9, rows with >= 1 lepton below 0.9, share with
//      kIsLepton == 0, exact polygon fractions (TCutG::IsInside on p_rec, beta_new from scan_nt) for the own species
//      at 1 sigma, 2 sigma and the PAT level (p, pi: 3 sigma, e: 5 sigma).
//   5. Protons in the pi+ slot with a DIFFERENT estimator than the chapter (which used the midpoint of the pion and
//      proton DST velocities): DST mass^2 = p_rec^2 (1/beta^2 - 1) > 0.35 (GeV/c^2)^2, pi- slot as the pion reference,
//      share = (f+ - f-) / (1 - f-) per momentum bin (500-1000-1500-1900-2200-2500 MeV/c, uncorrected p).
//   6. Writes the class of every FAT row (us, pp, mm): clsB = at least one lepton with DST beta < 0.9,
//      clsK = at least one lepton with kIsLepton == 0 (the sorter flag; the alternative definition), for
//      verify_pid_peak.C.
// Inputs : /hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list (PAT day files, read-only),
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/output_pippim{epem,epep,emem}_exp.root,
//          polygons: /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/remote/pid_cuts/pp45_pid_cuts_exp_betanew_v1.root
//          (the copy taken from vae25, md5 b8f0cc2f00c3c0e6d50627ea0bfc8c22, not the chapter's inputs/ copy)
// Outputs: /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/verify/verify_pid_rows.txt
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/verify/verify_pid_class.root
// Run:
//   mkdir -p /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/verify
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_pid_verify", kTRUE)' \
//        '/home/przygoda/HADES/FAT/FAT/plots/exp_report/verify_pid_rows.C+'
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TString.h>
#include <TSystem.h>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace vpr {
const char* kList  = "/hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list";
const char* kFat   = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/output_pippim%s_exp.root";
const char* kCuts  = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/remote/pid_cuts/pp45_pid_cuts_exp_betanew_v1.root";
const char* kOut   = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/verify/";

struct PatRow { float b[4]; float tof[2]; };      // DST beta of pip, pim, l1, l2; DST tof of l1, l2
}

void verify_pid_rows()
{
   using namespace vpr;
   gSystem->mkdir(kOut, kTRUE);
   FILE* fo = fopen(Form("%sverify_pid_rows.txt", kOut), "w");
   auto out = [&](const char* s) { printf("%s", s); fprintf(fo, "%s", s); };

   // ---- polygons
   TFile fc(kCuts, "READ");
   const char* sp[4] = {"pip", "pim", "ep", "em"};
   TCutG *c1[4], *c2[4], *cP[4];
   for (int k = 0; k < 4; ++k) {
      c1[k] = (TCutG*)fc.Get(Form("%s_cut_1sig", sp[k]));
      c2[k] = (TCutG*)fc.Get(Form("%s_cut_2sig", sp[k]));
      cP[k] = (TCutG*)fc.Get(Form("%s_cut_%s", sp[k], k < 2 ? "3sig" : "5sig"));
   }

   // ---- 1. DST beta from the PAT day files
   std::vector<std::string> files;
   { std::ifstream in(kList); std::string l; while (std::getline(in, l)) if (!l.empty() && l[0] != '#') files.push_back(l); }
   const char* tree[3] = {"PipPimEpEm_ID", "PipPimEpEp_ID", "PipPimEmEm_ID"};
   const char* l1[3] = {"ep", "ep1", "em1"};
   const char* l2[3] = {"em", "ep2", "em2"};
   std::map<int, std::vector<PatRow>> pat[3];
   for (auto& fn : files) {
      const int day = atoi(TString(gSystem->BaseName(fn.c_str()))(4, 3).Data());
      TFile* f = TFile::Open(fn.c_str(), "READ");
      for (int s = 0; s < 3; ++s) {
         TTree* t = (TTree*)f->Get(tree[s]);
         t->SetBranchStatus("*", 0);
         PatRow r;
         auto bind = [&](const char* n, float* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
         bind("pip_beta", &r.b[0]); bind("pim_beta", &r.b[1]);
         bind(Form("%s_beta", l1[s]), &r.b[2]); bind(Form("%s_beta", l2[s]), &r.b[3]);
         bind(Form("%s_tof", l1[s]), &r.tof[0]); bind(Form("%s_tof", l2[s]), &r.tof[1]);
         std::vector<PatRow>& v = pat[s][day];
         v.resize(t->GetEntries());
         for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); v[i] = r; }
      }
      f->Close();
      delete f;
      printf("read day %d\n", day);
   }

   // ---- 2. FAT rows
   TFile fcl(Form("%sverify_pid_class.root", kOut), "RECREATE");
   const char* fatn[3] = {"epem", "epep", "emem"};
   const char* cn[3] = {"cls_us", "cls_pp", "cls_mm"};
   const int NS = 8;
   const char* sn[NS] = {"fat", "pre", "pre_best", "pre_vtx", "peak", "bump", "peak_oa2", "bump_oa2"};
   for (int s = 0; s < 3; ++s) {
      TFile* f = TFile::Open(Form(kFat, fatn[s]), "READ");
      TTree* t = (TTree*)f->Get("scan_nt");
      t->SetBranchStatus("*", 0);
      float fid, pe, cor, vz, isb, uvp, uvm, oal, mpp, m4, oa, kl1, kl2;
      float prec[4], bn[4];
      auto bind = [&](const char* n, float* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
      bind("file_id", &fid); bind("pat_entry", &pe); bind("cor_ok", &cor); bind("eVertReco_z", &vz); bind("isBest", &isb);
      bind("pip_isUsedVertex", &uvp); bind("pim_isUsedVertex", &uvm); bind("oa_pippim_epem_lab_cor", &oal);
      bind("m_pippim_cor", &mpp); bind("m_pippimepem_cor", &m4); bind("oa", &oa);
      bind("ep_kIsLepton", &kl1); bind("em_kIsLepton", &kl2);
      for (int k = 0; k < 4; ++k) { bind(Form("%s_p_rec", sp[k]), &prec[k]); bind(Form("%s_beta_new", sp[k]), &bn[k]); }

      fcl.cd();                                           // the class tree lives in the output file
      TTree* tc = new TTree(cn[s], "class per FAT row");
      int oFid, oPe, oB, oK;
      tc->Branch("file_id", &oFid, "file_id/I"); tc->Branch("pat_entry", &oPe, "pat_entry/I");
      tc->Branch("clsB", &oB, "clsB/I"); tc->Branch("clsK", &oK, "clsK/I");

      double n[NS] = {0}, low1[NS] = {0}, low2[NS] = {0}, low12[NS] = {0}, lowBoth[NS] = {0}, kl0a[NS] = {0}, kl0b[NS] = {0};
      std::map<int, double> dayN, dayLow1, dayLow2;
      long long tof01 = 0, tof1 = 0, tofNeg = 0, tof60 = 0, tof4060 = 0;
      double in1[NS][4] = {{0}}, in2[NS][4] = {{0}}, inP[NS][4] = {{0}};
      // protons: bins x (pi+ all, pi+ heavy, pi- all, pi- heavy) for sel fat, pre, peak, bump
      const double edges[6] = {500, 1000, 1500, 1900, 2200, 2500};
      double pr[NS][5][4] = {{{0}}};
      long long nMiss = 0, nZero = 0;
      double scaleM = 0;
      const Long64_t nent = t->GetEntries();
      for (Long64_t i = 0; i < nent; ++i) {
         t->GetEntry(i);
         if (i == 0) scaleM = (m4 > 10.f) ? 1000. : 1.;     // MeV or GeV
         const int day = (int)std::lround(fid), ent = (int)std::lround(pe);
         auto it = pat[s].find(day);
         if (it == pat[s].end() || ent < 0 || ent >= (int)it->second.size()) { ++nMiss; continue; }
         const PatRow& r = it->second[ent];
         const bool B = r.b[2] < 0.9f || r.b[3] < 0.9f;
         const bool K = std::lround(kl1) == 0 || std::lround(kl2) == 0;
         oFid = day; oPe = ent; oB = B; oK = K; tc->Fill();
         if (s != 0) { n[0] += 1; if (B) low12[0] += 1; continue; }
         const double M4 = m4 / scaleM, Mpp = mpp / scaleM;
         bool sel[NS];
         sel[0] = true;
         sel[1] = std::lround(cor) == 1 && vz > -160.f;
         sel[2] = sel[1] && std::lround(isb) == 1;
         sel[3] = sel[1] && std::lround(uvp) == 1 && std::lround(uvm) == 1;
         sel[4] = sel[1] && oal < 55.f && Mpp > 0.30 && M4 >= 0.53364 && M4 < 0.55416;
         sel[5] = sel[1] && oal < 40.f && Mpp > 0.30 && Mpp < 0.41 && M4 >= 0.4550 && M4 < 0.5225;
         sel[6] = sel[4] && oa > 2.f;
         sel[7] = sel[5] && oa > 2.f;
         bool i1[4], i2[4], iP[4];
         for (int k = 0; k < 4; ++k) {
            i1[k] = c1[k]->IsInside(prec[k], bn[k]);
            i2[k] = c2[k]->IsInside(prec[k], bn[k]);
            iP[k] = cP[k]->IsInside(prec[k], bn[k]);
         }
         if (r.b[2] == 0.f && r.b[3] == 0.f) ++nZero;
         // per day (fat rows) and lepton DST tof near zero
         dayN[day] += 1; dayLow1[day] += r.b[2] < 0.9f; dayLow2[day] += r.b[3] < 0.9f;
         for (int j = 0; j < 2; ++j) {
            if (r.tof[j] >= 0.f && r.tof[j] < 0.1f) ++tof01;
            if (r.tof[j] >= 0.f && r.tof[j] < 1.f) ++tof1;
            if (r.tof[j] < 0.f) ++tofNeg;
            if (r.tof[j] > 60.f) ++tof60;
            if (r.tof[j] >= 40.f && r.tof[j] < 60.f) ++tof4060;
         }
         for (int q = 0; q < NS; ++q) {
            if (!sel[q]) continue;
            n[q] += 1;
            if (r.b[2] < 0.9f) low1[q] += 1;
            if (r.b[3] < 0.9f) low2[q] += 1;
            if (B) low12[q] += 1;
            if (r.b[2] < 0.9f && r.b[3] < 0.9f) lowBoth[q] += 1;
            if (std::lround(kl1) == 0) kl0a[q] += 1;
            if (std::lround(kl2) == 0) kl0b[q] += 1;
            for (int k = 0; k < 4; ++k) { in1[q][k] += i1[k]; in2[q][k] += i2[k]; inP[q][k] += iP[k]; }
            for (int k = 0; k < 2; ++k) {
               const double p = prec[k], b = r.b[k];
               for (int j = 0; j < 5; ++j) {
                  if (!(p >= edges[j] && p < edges[j + 1])) continue;
                  pr[q][j][2 * k] += 1;
                  const double m2 = (b > 0) ? (p / 1000.) * (p / 1000.) * (1. / (b * b) - 1.) : 99.;
                  if (m2 > 0.35) pr[q][j][2 * k + 1] += 1;
               }
            }
         }
      }
      out(Form("\n=== %s (FAT scan_nt %s): %lld rows, not matched to PAT %lld, rows with both lepton DST beta == 0: %lld, mass unit scale %.0f\n",
               cn[s], fatn[s], (long long)nent, nMiss, nZero, scaleM));
      if (s != 0) {
         out(Form("  like-sign rows %.0f: rows with >= 1 lepton DST beta < 0.9: %.4f\n", n[0], low12[0] / n[0]));
      } else {
         for (int q = 0; q < NS; ++q) {
            const double N = n[q];
            out(Form("  %-9s rows %8.0f | e+ beta<0.9 %.4f +- %.4f | e- %.4f +- %.4f | rows >=1 below %.4f, both %.4f | kIsLepton==0 e+ %.4f e- %.4f\n",
                     sn[q], N, low1[q] / N, std::sqrt(low1[q] / N * (1 - low1[q] / N) / N), low2[q] / N, std::sqrt(low2[q] / N * (1 - low2[q] / N) / N),
                     low12[q] / N, lowBoth[q] / N, kl0a[q] / N, kl0b[q] / N));
            out(Form("            inside own 1sig pi+ %.4f pi- %.4f e+ %.4f e- %.4f | 2sig %.4f %.4f %.4f %.4f | PAT level %.6f %.6f %.6f %.6f\n",
                     in1[q][0] / N, in1[q][1] / N, in1[q][2] / N, in1[q][3] / N, in2[q][0] / N, in2[q][1] / N, in2[q][2] / N, in2[q][3] / N,
                     inP[q][0] / N, inP[q][1] / N, inP[q][2] / N, inP[q][3] / N));
            if (q == 0 || q == 1 || q == 4 || q == 5) {
               double tot = 0;
               for (int j = 0; j < 5; ++j) {
                  const double fp = pr[q][j][0] > 0 ? pr[q][j][1] / pr[q][j][0] : 0, fm = pr[q][j][2] > 0 ? pr[q][j][3] / pr[q][j][2] : 0;
                  const double sh = pr[q][j][2] > 20 ? (fp - fm) / (1 - fm) : fp;
                  tot += sh * pr[q][j][0];
                  out(Form("            protons in pi+ (DST m2 > 0.35): p %4.0f-%4.0f  pi+ %8.0f f+ %.4f  pi- %8.0f f- %.4f  share %.3f  protons %.0f\n",
                           edges[j], edges[j + 1], pr[q][j][0], fp, pr[q][j][2], fm, sh, sh * pr[q][j][0]));
               }
               out(Form("            protons in pi+ total %.0f = %.4f of the rows\n", tot, tot / N));
            }
         }
      }
      if (s == 0) {
         out(Form("  lepton DST tof (fat rows, both slots): [0,0.1) ns %lld, [0,1) ns %lld, < 0 %lld, > 60 ns %lld, 40-60 ns %lld\n", tof01, tof1, tofNeg, tof60, tof4060));
         // periods: I 034-046, IIA 047-048, IIB 049-053, III 059-067
         const int plo[4] = {34, 47, 49, 59}, phi[4] = {46, 48, 53, 67};
         const char* pn[4] = {"I", "IIA", "IIB", "III"};
         for (int q = 0; q < 4; ++q) {
            double a = 0, b1 = 0, b2 = 0;
            for (auto& kv : dayN) if (kv.first >= plo[q] && kv.first <= phi[q]) { a += kv.second; b1 += dayLow1[kv.first]; b2 += dayLow2[kv.first]; }
            out(Form("  period %-3s rows %8.0f  e+ beta<0.9 %.4f  e- %.4f\n", pn[q], a, b1 / a, b2 / a));
         }
         int dmin = -1, dmax = -1; double smin = 9, smax = -1;
         for (auto& kv : dayN) {
            const double sh = 0.5 * (dayLow1[kv.first] + dayLow2[kv.first]) / kv.second;
            if (sh < smin) { smin = sh; dmin = kv.first; }
            if (sh > smax) { smax = sh; dmax = kv.first; }
         }
         out(Form("  per day (mean of e+ and e-): min %.4f (day %03d: e+ %.4f e- %.4f), max %.4f (day %03d: e+ %.4f e- %.4f)\n", smin, dmin,
                  dayLow1[dmin] / dayN[dmin], dayLow2[dmin] / dayN[dmin], smax, dmax, dayLow1[dmax] / dayN[dmax], dayLow2[dmax] / dayN[dmax]));
      }
      fcl.cd();
      tc->Write();
      f->Close();
      delete f;
   }
   fcl.Close();
   fclose(fo);
}
