// pid_final_fill.C -- PID chapter of the experimental-data report: tracks of the FINAL analysed sample.
//
// What it does
//   Reads ONE day file of the PAT betanew production (read-only), tree PipPimEpEm_ID (and the like-sign
//   trees PipPimEpEp_ID, PipPimEmEm_ID for the tag trees), keeps the rows that FAT keeps
//   (isBest >= 0, eVertReco_z > -500 mm, trigbit == 8192 = PT3) and, for the four track slots
//   pip (+), pim (-), ep (+), em (-), fills
//     - beta_new versus charge x momentum (x = sign x <slot>_p, the uncorrected momentum used by the
//       PID polygons; same binning as the simulation report: 10 MeV/c x 0.0025),
//     - the DST beta (START based, HParticleCand::getBeta) versus charge x momentum,
//     - for the lepton slots: 1D DST beta (split by the sorter flag kIsLepton), the time offset
//       dt = tof - L/c (tof = DST time of flight against the START time, L = track length), per-day counts,
//     - counters: track entries (row-slots) and unique tracks (each track once per event group:
//       adjacent rows with equal (runnumber, event) inside the day file), the fraction inside the
//       betanew_v1 polygon of the slot's own species at 1, 2, 2.5, 3, 3.5, 5 sigma (exact
//       TCutG::IsInside on the track values) and inside the PAT-level polygons of the other species
//       of the same charge (p, pi: 3 sigma; e: 5 sigma), DST beta < 0.9, kIsLepton, evtStart.
//   Selections (rows), defined from the compact scan file (tree us, key file_id = day number, pat_entry):
//     0 fat      : FAT sample = isBest >= 0 && vz > -500 && PT3                  (all rows of the compact)
//     1 pre      : fat && cor_ok == 1 && vz > -160
//     2 pre_best : pre && isBest == 1
//     3 pre_vtx  : pre && pip_isUsedVertex && pim_isUsedVertex (both pions in the event-vertex fit)
//     4 peak     : pre && oa_lab < 55 && M(pipi) > 0.30 && 0.53364 <= M4 < 0.55416 (data peak window)
//     5 bump     : pre && oa_lab < 40 && 0.30 < M(pipi) < 0.41 && 0.4550 <= M4 < 0.5225 (data bump window)
//   (M4, M(pipi) from the energy-loss corrected momenta, as in plots/ScanEngine.h; oa_lab in deg.)
//   Tag trees tag_us / tag_pp / tag_mm: one entry per FAT row of each hypothesis with file_id, pat_entry,
//   the DST beta, tof, track length and kIsLepton of the two lepton slots (us: ep, em; pp: ep1, ep2;
//   mm: em1, em2) and evtStart; used by pid_peak_dstbeta.C to split the eta signal by lepton DST beta.
//
// Inputs
//   patFile : PAT day file, e.g. /hdd3/przygoda/HADES/files/exp/gen4/lep_049.root (from lep_betanew.list)
//   compact : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root
//   cutFile : PID polygons of the data production, pp45_pid_cuts_exp_betanew_v1.root
//   out     : output ROOT file (one per day; hadd them afterwards)
//   maxrows : read at most this many rows of PipPimEpEm_ID (0 = all; for tests)
// Output
//   ROOT file with the histograms above (names hBN_<sel>_<slot>, hBD_<sel>_<slot>, hDst*_<sel>_<slot>,
//   hDt_<sel>_<slot>, hDay*_<sel>_<slot>, hCnt_<sel>_<slot>, hRows) and the tag trees.
//
// Run (one day; the build directory is outside the repository):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_pid", kTRUE)' \
//     '/home/przygoda/HADES/FAT/FAT/plots/exp_report/pid_final_fill.C+("/hdd3/przygoda/HADES/files/exp/gen4/lep_049.root", "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root", "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/inputs/pp45_pid_cuts_exp_betanew_v1.root", "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/pid/final/pid_final_049.root")'
//   All 29 days (6 in parallel) + hadd: plots/exp_report/pid_run_final.sh
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TString.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <tuple>
#include <vector>

namespace pidff {
const int kNSel = 6;
const char* kSel[kNSel] = {"fat", "pre", "pre_best", "pre_vtx", "peak", "bump"};
const int kNSlot = 4;
const char* kSlot[kNSlot] = {"pip", "pim", "ep", "em"};
const int kSign[kNSlot] = {+1, -1, +1, -1};
const int kNLev = 6;
const char* kLev[kNLev] = {"1sig", "2sig", "25sig", "3sig", "35sig", "5sig"};
// counters (bin labels of hCnt_<sel>_<slot>)
enum {
   cN = 1, cNU,                                          // track entries, unique tracks
   cIn1, cIn2, cIn25, cIn3, cIn35, cIn5,                 // inside own polygon at the 6 levels (entries)
   cInU1, cInU2, cInU25, cInU3, cInU35, cInU5,           // the same for unique tracks
   cInP, cInPip, cInPim, cInEp, cInEm,                   // inside the PAT-level polygon of p, pi+, pi-, e+, e- (entries)
   cLow, cLowU,                                          // DST beta < 0.9 (entries, unique)
   cKL0, cLowKL0, cLowKL1,                               // kIsLepton == 0; beta<0.9 && kIsLepton==0; beta<0.9 && kIsLepton==1
   cGt1, cStart0, cLowStart0, cRing, cBneg, cTof60,      // DST beta > 1; evtStart == 0; beta<0.9 && evtStart==0; isring != 0; beta < 0; tof > 60
   cPgt19, cNBin
};
const char* kCntLab[cNBin] = {"", "n", "n_unique", "in_1sig", "in_2sig", "in_25sig", "in_3sig", "in_35sig", "in_5sig",
                              "inU_1sig", "inU_2sig", "inU_25sig", "inU_3sig", "inU_35sig", "inU_5sig",
                              "in_p3", "in_pip3", "in_pim3", "in_ep5", "in_em5",
                              "dstbeta_lt09", "dstbeta_lt09_unique", "kIsLepton0", "lt09_kIsLepton0", "lt09_kIsLepton1",
                              "dstbeta_gt1", "evtStart0", "lt09_evtStart0", "isring", "dstbeta_lt0", "tof_gt60", "p_gt1900"};
}

void pid_final_fill(const char* patFile, const char* compact, const char* cutFile, const char* out, long long maxrows = 0)
{
   using namespace pidff;
   TStopwatch sw;
   // ---- day number from the file name lep_XXX.root
   TString base = gSystem->BaseName(patFile);
   const int day = atoi(TString(base(4, 3)).Data());
   printf("pid_final_fill: %s day %d\n", patFile, day);

   // ---- polygons
   TFile fc(cutFile, "READ");
   const char* sp[5] = {"p", "pip", "pim", "ep", "em"};
   TCutG* cut[5][kNLev];
   for (int s = 0; s < 5; ++s)
      for (int l = 0; l < kNLev; ++l) {
         cut[s][l] = (TCutG*)fc.Get(Form("%s_cut_%s", sp[s], kLev[l]));
         if (!cut[s][l]) { printf("missing polygon %s_cut_%s\n", sp[s], kLev[l]); return; }
      }
   const int kOwn[kNSlot] = {1, 2, 3, 4};       // species index of the slot
   const int kPat[5] = {3, 3, 3, 5, 5};         // PAT level index (3sig = 3, 5sig = 5) per species

   // ---- selection flags from the compact file (rows of this day)
   std::vector<unsigned char> flag;            // bit 0 present, 1..5 = selections 1..5
   long long nCompact = 0;
   {
      TFile f(compact, "READ");
      TTree* t = (TTree*)f.Get("us");
      int fid, pe, cor, isb, uvp, uvm;
      float vz, oal, mpp, m4;
      t->SetBranchStatus("*", 0);
      for (const char* b : {"file_id", "pat_entry", "cor_ok", "isBest", "pip_isUsedVertex", "pim_isUsedVertex", "vz", "oa_lab", "mpipi", "m4"})
         t->SetBranchStatus(b, 1);
      t->SetBranchAddress("file_id", &fid); t->SetBranchAddress("pat_entry", &pe); t->SetBranchAddress("cor_ok", &cor);
      t->SetBranchAddress("isBest", &isb); t->SetBranchAddress("pip_isUsedVertex", &uvp); t->SetBranchAddress("pim_isUsedVertex", &uvm);
      t->SetBranchAddress("vz", &vz); t->SetBranchAddress("oa_lab", &oal); t->SetBranchAddress("mpipi", &mpp); t->SetBranchAddress("m4", &m4);
      const Long64_t n = t->GetEntries();
      for (Long64_t i = 0; i < n; ++i) {
         t->GetEntry(i);
         if (fid != day) continue;
         if (pe < 0) continue;
         if ((size_t)pe >= flag.size()) flag.resize(pe + 1000, 0);
         unsigned char b = 1;
         const bool pre = cor == 1 && vz > -160.f;
         if (pre) b |= 1 << 1;
         if (pre && isb == 1) b |= 1 << 2;
         if (pre && uvp == 1 && uvm == 1) b |= 1 << 3;
         if (pre && oal < 55.f && mpp > 0.30f && m4 >= 0.53364f && m4 < 0.55416f) b |= 1 << 4;
         if (pre && oal < 40.f && mpp > 0.30f && mpp < 0.41f && m4 >= 0.4550f && m4 < 0.5225f) b |= 1 << 5;
         flag[pe] = b;
         ++nCompact;
      }
      f.Close();
   }
   printf("  compact rows of day %d: %lld (%.1f s)\n", day, nCompact, sw.RealTime()); sw.Continue();

   // ---- output histograms
   TFile fo(out, "RECREATE");
   const int nx = 600, ny = 520; const double x0 = -2500, x1 = 3500;
   TH2D *hBN[kNSel][kNSlot], *hBD[kNSel][kNSlot];
   TH1D *hCnt[kNSel][kNSlot], *hDst[kNSel][kNSlot], *hDstL[kNSel][kNSlot], *hDstH[kNSel][kNSlot], *hDt[kNSel][kNSlot];
   TH1D *hDay[kNSel][kNSlot], *hDayLow[kNSel][kNSlot];
   for (int s = 0; s < kNSel; ++s)
      for (int k = 0; k < kNSlot; ++k) {
         const char* sn = kSel[s]; const char* kn = kSlot[k];
         hBN[s][k] = new TH2D(Form("hBN_%s_%s", sn, kn), Form("%s %s;q p [MeV/c];#beta_{new}", sn, kn), nx, x0, x1, ny, 0., 1.3);
         hBD[s][k] = new TH2D(Form("hBD_%s_%s", sn, kn), Form("%s %s;q p [MeV/c];#beta (DST, START)", sn, kn), nx, x0, x1, ny, 0., 1.3);
         hCnt[s][k] = new TH1D(Form("hCnt_%s_%s", sn, kn), Form("counters %s %s", sn, kn), cNBin - 1, 0.5, cNBin - 0.5);
         for (int b = 1; b < cNBin; ++b) hCnt[s][k]->GetXaxis()->SetBinLabel(b, kCntLab[b]);
         hDst[s][k] = new TH1D(Form("hDst_%s_%s", sn, kn), Form("%s %s;#beta (DST, START);tracks", sn, kn), 340, -0.1, 1.6);
         hDstL[s][k] = new TH1D(Form("hDstL_%s_%s", sn, kn), Form("%s %s kIsLepton=1;#beta (DST, START);tracks", sn, kn), 340, -0.1, 1.6);
         hDstH[s][k] = new TH1D(Form("hDstH_%s_%s", sn, kn), Form("%s %s kIsLepton=0;#beta (DST, START);tracks", sn, kn), 340, -0.1, 1.6);
         hDt[s][k] = new TH1D(Form("hDt_%s_%s", sn, kn), Form("%s %s;t_{DST} - L/c [ns];tracks", sn, kn), 900, -20., 70.);
         hDay[s][k] = new TH1D(Form("hDay_%s_%s", sn, kn), Form("%s %s;day;tracks", sn, kn), 40, 33.5, 73.5);
         hDayLow[s][k] = new TH1D(Form("hDayLow_%s_%s", sn, kn), Form("%s %s DST beta<0.9;day;tracks", sn, kn), 40, 33.5, 73.5);
      }
   TH1D* hRows = new TH1D("hRows", "rows", 8, 0.5, 8.5);
   const char* rl[8] = {"pat_rows", "pat_fat_rows", "compact_rows", "fat_rows_in_compact", "fat_rows_not_in_compact", "pp_fat_rows", "mm_fat_rows", "files"};
   for (int b = 0; b < 8; ++b) hRows->GetXaxis()->SetBinLabel(b + 1, rl[b]);
   hRows->Fill(3, (double)nCompact);
   hRows->Fill(8);

   // ---- tag trees
   int tFid = day, tPe = 0, tK1 = 0, tK2 = 0, tStart = 0;
   float tB1 = 0, tB2 = 0, tT1 = 0, tT2 = 0, tL1 = 0, tL2 = 0;
   TTree* tag[3];
   const char* tagName[3] = {"tag_us", "tag_pp", "tag_mm"};
   for (int i = 0; i < 3; ++i) {
      tag[i] = new TTree(tagName[i], "lepton DST beta of the FAT rows (isBest>=0, vz>-500, PT3)");
      tag[i]->Branch("file_id", &tFid, "file_id/I"); tag[i]->Branch("pat_entry", &tPe, "pat_entry/I");
      tag[i]->Branch("b1", &tB1, "b1/F"); tag[i]->Branch("b2", &tB2, "b2/F");
      tag[i]->Branch("tof1", &tT1, "tof1/F"); tag[i]->Branch("tof2", &tT2, "tof2/F");
      tag[i]->Branch("len1", &tL1, "len1/F"); tag[i]->Branch("len2", &tL2, "len2/F");
      tag[i]->Branch("k1", &tK1, "k1/I"); tag[i]->Branch("k2", &tK2, "k2/I"); tag[i]->Branch("evtStart", &tStart, "evtStart/I");
   }

   TFile* fin = TFile::Open(patFile, "READ");
   if (!fin || fin->IsZombie()) { printf("cannot open %s\n", patFile); return; }
   const double cLight = 299.792458;   // mm/ns

   // ---- like-sign trees: tags only
   {
      const char* tn[2] = {"PipPimEpEp_ID", "PipPimEmEm_ID"};
      const char* pf[2][2] = {{"ep1", "ep2"}, {"em1", "em2"}};
      for (int q = 0; q < 2; ++q) {
         TTree* t = (TTree*)fin->Get(tn[q]);
         if (!t) { printf("no %s\n", tn[q]); continue; }
         t->SetBranchStatus("*", 0);
         float isb, vz, trig, est, b[2], tof[2], len[2], kl[2];
         auto use = [&](const char* n, float* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
         use("isBest", &isb); use("eVertReco_z", &vz); use("trigbit", &trig); use("evtStart", &est);
         for (int j = 0; j < 2; ++j) {
            use(Form("%s_beta", pf[q][j]), &b[j]); use(Form("%s_tof", pf[q][j]), &tof[j]);
            use(Form("%s_tracklength", pf[q][j]), &len[j]); use(Form("%s_kIsLepton", pf[q][j]), &kl[j]);
         }
         const Long64_t n = t->GetEntries();
         long long nf = 0;
         for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            if (!(isb > -0.5f && vz > -500.f && std::lround(trig) == 8192)) continue;
            tPe = (int)i; tB1 = b[0]; tB2 = b[1]; tT1 = tof[0]; tT2 = tof[1]; tL1 = len[0]; tL2 = len[1];
            tK1 = (int)std::lround(kl[0]); tK2 = (int)std::lround(kl[1]); tStart = (int)std::lround(est);
            tag[q + 1]->Fill();
            ++nf;
         }
         hRows->Fill(6 + q, (double)nf);
         printf("  %s: %lld rows, FAT rows %lld\n", tn[q], (long long)n, nf);
      }
   }

   // ---- unlike-sign tree
   TTree* t = (TTree*)fin->Get("PipPimEpEm_ID");
   t->SetBranchStatus("*", 0);
   float isb, vz, trig, est, run, evt;
   float p[kNSlot], bn[kNSlot], bd[kNSlot], tof[kNSlot], len[kNSlot], kl[kNSlot], ring[kNSlot], th[kNSlot], ph[kNSlot];
   auto use = [&](const char* n, float* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
   use("isBest", &isb); use("eVertReco_z", &vz); use("trigbit", &trig); use("evtStart", &est);
   use("runnumber", &run); use("event", &evt);
   for (int k = 0; k < kNSlot; ++k) {
      use(Form("%s_p", kSlot[k]), &p[k]); use(Form("%s_beta_new", kSlot[k]), &bn[k]); use(Form("%s_beta", kSlot[k]), &bd[k]);
      use(Form("%s_tof", kSlot[k]), &tof[k]); use(Form("%s_tracklength", kSlot[k]), &len[k]);
      use(Form("%s_kIsLepton", kSlot[k]), &kl[k]); use(Form("%s_isring", kSlot[k]), &ring[k]);
      use(Form("%s_theta", kSlot[k]), &th[k]); use(Form("%s_phi", kSlot[k]), &ph[k]);
   }
   Long64_t n = t->GetEntries();
   if (maxrows > 0 && n > maxrows) n = maxrows;
   hRows->Fill(1, (double)n);
   float lastRun = -1e9, lastEvt = -1e9;
   std::set<std::tuple<int, int, long, long, long>> seen[kNSel];
   long long nFat = 0, nIn = 0, nOut = 0;
   for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      if (run != lastRun || evt != lastEvt) { for (int s = 0; s < kNSel; ++s) seen[s].clear(); lastRun = run; lastEvt = evt; }
      if (!(isb > -0.5f && vz > -500.f && std::lround(trig) == 8192)) continue;
      ++nFat;
      // tag of the unlike-sign row
      tPe = (int)i; tB1 = bd[2]; tB2 = bd[3]; tT1 = tof[2]; tT2 = tof[3]; tL1 = len[2]; tL2 = len[3];
      tK1 = (int)std::lround(kl[2]); tK2 = (int)std::lround(kl[3]); tStart = (int)std::lround(est);
      tag[0]->Fill();
      const unsigned char fl = (i < (Long64_t)flag.size()) ? flag[i] : 0;
      if (fl & 1) ++nIn; else ++nOut;
      bool sel[kNSel];
      sel[0] = true;
      for (int s = 1; s < kNSel; ++s) sel[s] = (fl >> s) & 1;
      for (int k = 0; k < kNSlot; ++k) {
         const double x = kSign[k] * p[k];
         // polygon tests once per track
         bool inOwn[kNLev];
         const int own = kOwn[k];
         for (int l = 0; l < kNLev; ++l) inOwn[l] = cut[own][l]->IsInside(p[k], bn[k]);
         bool inSp[5] = {false, false, false, false, false};
         for (int sp5 = 0; sp5 < 5; ++sp5) {
            const bool pos = (sp5 == 0 || sp5 == 1 || sp5 == 3);
            if ((kSign[k] > 0) != pos) continue;
            inSp[sp5] = (sp5 == own) ? inOwn[kPat[sp5]] : cut[sp5][kPat[sp5]]->IsInside(p[k], bn[k]);
         }
         const bool low = bd[k] < 0.9f;
         const bool kl0 = std::lround(kl[k]) == 0;
         const bool st0 = std::lround(est) == 0;
         const auto key = std::make_tuple(k, kSign[k], std::lround(p[k] * 10.), std::lround(th[k] * 1000.), std::lround(ph[k] * 1000.));
         const double dt = tof[k] - len[k] / cLight;
         for (int s = 0; s < kNSel; ++s) {
            if (!sel[s]) continue;
            TH1D* c = hCnt[s][k];
            const bool uniq = seen[s].insert(key).second;
            hBN[s][k]->Fill(x, bn[k]);
            hBD[s][k]->Fill(x, bd[k]);
            c->Fill(cN);
            if (uniq) c->Fill(cNU);
            for (int l = 0; l < kNLev; ++l) if (inOwn[l]) { c->Fill(cIn1 + l); if (uniq) c->Fill(cInU1 + l); }
            for (int sp5 = 0; sp5 < 5; ++sp5) if (inSp[sp5]) c->Fill(cInP + sp5);
            if (low) { c->Fill(cLow); if (uniq) c->Fill(cLowU); }
            if (kl0) c->Fill(cKL0);
            if (low && kl0) c->Fill(cLowKL0);
            if (low && !kl0) c->Fill(cLowKL1);
            if (bd[k] > 1.f) c->Fill(cGt1);
            if (st0) c->Fill(cStart0);
            if (low && st0) c->Fill(cLowStart0);
            if (ring[k] > 0.5f) c->Fill(cRing);
            if (bd[k] < 0.f) c->Fill(cBneg);
            if (tof[k] > 60.f) c->Fill(cTof60);
            if (p[k] > 1900.f) c->Fill(cPgt19);
            hDst[s][k]->Fill(bd[k]);
            (kl0 ? hDstH : hDstL)[s][k]->Fill(bd[k]);
            hDt[s][k]->Fill(dt);
            hDay[s][k]->Fill(day);
            if (low) hDayLow[s][k]->Fill(day);
         }
      }
   }
   hRows->Fill(2, (double)nFat);
   hRows->Fill(4, (double)nIn);
   hRows->Fill(5, (double)nOut);
   printf("  PipPimEpEm_ID: %lld rows read, FAT rows %lld, of them in compact %lld, not in compact %lld (compact rows %lld) %.1f s\n",
          (long long)n, nFat, nIn, nOut, nCompact, sw.RealTime());
   fin->Close();
   fo.cd();
   fo.Write();
   fo.Close();
   fc.Close();
}
