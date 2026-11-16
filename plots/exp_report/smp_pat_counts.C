// ============================================================================
// smp_pat_counts.C  --  row, event and run counts of ONE day file of the PAT
//                       input (first steps of the cut flow), data
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Input: /hdd3/przygoda/HADES/files/exp/gen4/lep_<day>.root (PAT production tag
// prod-feb22-betanew-2026-10-28, Hydra 2.7.6, PID polygons betanew_v1), the file of
// the list lep_betanew.list that FAT reads. Only the entry counts and a few
// branches are read (the files are 5-21 GB):
//   PipPimEpEm_ID, PipPimEpEp_ID, PipPimEmEm_ID: isBest, eVertReco_z, trigbit,
//                                               runnumber, event
//   EpEm_ID: runnumber only (runs present in the day file, as in the luminosity
//            coverage check)
// What it does, per hypothesis tree (one row = one stored combination):
//   steps of the FAT event selection, in the order of main.cc processEvent:
//     s0 all rows; s1 isBest >= 0 (FAT "isBestAll"); s2 + eVertReco_z > -500 mm;
//     s3 + trigbit == 8192 (PT3)  ->  s3 = the rows FAT writes
//   per step: rows, rows with isBest == 1, event groups (maximal runs of adjacent
//   rows with equal (runnumber, event) inside the day file) with at least one row
//   passing; trigbit values before and after the vertex cut; histograms of
//   eVertReco_z (all rows, and isBest == 1 rows = one per event); rows per run.
// Outputs (directory /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/smp/pat):
//   pat_<day>.txt   key-value lines (read by smp_report.C)
//   pat_<day>.root  vz histograms
//   runs_<day>.txt  run, EpEm_ID rows, rows of the three trees at s0 and s3
//
// Run (from the repository root /home/przygoda/HADES/FAT/FAT), one day, e.g. 53:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_smp", kTRUE)' 'plots/exp_report/smp_pat_counts.C+(53)'
// All 29 days (6 in parallel; compile once first with the line above):
//   for d in 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 59 60 61 62 63 64 65 66 67; do echo $d; done | \
//   xargs -P 6 -I{} root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_smp", kTRUE)' 'plots/exp_report/smp_pat_counts.C+({})'
// ============================================================================
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"

#include "smp_common.h"

void smp_pat_counts(int day)
{
    TStopwatch sw;
    const std::string path = Form("%s/lep_%03d.root", smp::kPatDir, day);
    const std::string outDir = std::string(smp::kWork) + "/pat";
    gSystem->mkdir(outDir.c_str(), kTRUE);

    FileStat_t fs;
    Long64_t fsize = -1;
    if (gSystem->GetPathInfo(path.c_str(), fs) == 0) fsize = fs.fSize;

    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) { printf("cannot open %s\n", path.c_str()); return; }

    FILE* o = fopen(Form("%s/pat_%03d.txt", outDir.c_str(), day), "w");
    fprintf(o, "# smp_pat_counts.C: PAT counts of %s\n", path.c_str());
    fprintf(o, "FILE %d %lld %s\n", day, fsize, path.c_str());

    TFile* fh = TFile::Open(Form("%s/pat_%03d.root", outDir.c_str(), day), "RECREATE");

    // run -> rows: [0] EpEm_ID, then per hypothesis (s0, s3)
    std::map<int, std::vector<Long64_t>> runRows;
    auto runEntry = [&](int run) -> std::vector<Long64_t>& {
        auto it = runRows.find(run);
        if (it == runRows.end()) it = runRows.insert(std::make_pair(run, std::vector<Long64_t>(7, 0))).first;
        return it->second;
    };

    // ---- EpEm_ID: runs present ----
    {
        TTree* t = (TTree*)f->Get("EpEm_ID");
        Long64_t n = t ? t->GetEntries() : -1;
        if (t) {
            t->SetBranchStatus("*", 0);
            t->SetBranchStatus("runnumber", 1);
            Float_t run = 0;
            t->SetBranchAddress("runnumber", &run);
            t->SetCacheSize(100 * 1024 * 1024);
            t->AddBranchToCache("runnumber", kTRUE);
            float last = -1; int lastRun = -1;
            Long64_t cnt = 0;
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                if (run != last) {
                    if (lastRun >= 0) runEntry(lastRun)[0] += cnt;
                    cnt = 0; last = run; lastRun = (int)std::lround(run);
                }
                cnt++;
            }
            if (lastRun >= 0) runEntry(lastRun)[0] += cnt;
            t->ResetBranchAddresses();
        }
        fprintf(o, "TREE EpEm_ID %lld\n", n);
    }

    // ---- the three pi+pi-l+l- hypotheses ----
    const char* trees[3] = {"PipPimEpEm_ID", "PipPimEpEp_ID", "PipPimEmEm_ID"};
    const char* tag[3] = {"us", "pp", "mm"};
    for (int k = 0; k < 3; ++k) {
        TTree* t = (TTree*)f->Get(trees[k]);
        if (!t) { fprintf(o, "TREE %s -1\n", trees[k]); continue; }
        const Long64_t n = t->GetEntries();
        fprintf(o, "TREE %s %lld\n", trees[k], n);
        t->SetBranchStatus("*", 0);
        Float_t isBest = 0, vz = 0, trig = 0, run = 0, ev = 0;
        const char* br[5] = {"isBest", "eVertReco_z", "trigbit", "runnumber", "event"};
        Float_t* ad[5] = {&isBest, &vz, &trig, &run, &ev};
        for (int b = 0; b < 5; ++b) { t->SetBranchStatus(br[b], 1); t->SetBranchAddress(br[b], ad[b]); }
        t->SetCacheSize(100 * 1024 * 1024);
        for (int b = 0; b < 5; ++b) t->AddBranchToCache(br[b], kTRUE);

        fh->cd();
        TH1D* hVz = new TH1D(Form("vz_%s", tag[k]), "eVertReco_z, all rows;v_{z} [mm];rows", 1600, -1100, 500);
        TH1D* hVzB = new TH1D(Form("vzbest_%s", tag[k]), "eVertReco_z, isBest == 1 rows;v_{z} [mm];rows", 1600, -1100, 500);

        const int NS = 4;
        Long64_t rows[NS] = {0, 0, 0, 0}, best[NS] = {0, 0, 0, 0}, evts[NS] = {0, 0, 0, 0};
        bool gp[NS] = {false, false, false, false};
        Long64_t isBestNeg = 0, groupsBestNot1 = 0, vzNoVertex = 0, vzLow = 0;
        int nBestGrp = 0;
        std::map<int, Long64_t> trigBefore, trigAfter;   // trigbit -> rows (s1, s2)
        float lr = -1, le = -1;
        bool have = false;
        auto flush = [&]() {
            if (!have) return;
            for (int s = 0; s < NS; ++s) if (gp[s]) evts[s]++;
            if (nBestGrp != 1) groupsBestNot1++;
            for (int s = 0; s < NS; ++s) gp[s] = false;
            nBestGrp = 0;
        };
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            if (!have || run != lr || ev != le) { flush(); lr = run; le = ev; have = true; }
            const bool b1 = (isBest == 1.f);
            if (b1) nBestGrp++;
            std::vector<Long64_t>& rr = runEntry((int)std::lround(run));
            rr[1 + 2 * k] += 1;
            bool pass[NS];
            pass[0] = true;
            pass[1] = isBest > -0.5f;                       // MinCut strict, -0.5
            pass[2] = pass[1] && vz > -500.f;               // MinCut strict
            pass[3] = pass[2] && trig == 8192.f;            // ValueCut
            if (isBest < -0.5f) isBestNeg++;
            if (pass[1]) {
                hVz->Fill(vz);
                if (b1) hVzB->Fill(vz);
                if (vz < -999.f) vzNoVertex++;
                if (!(vz > -500.f)) vzLow++;
                trigBefore[(int)std::lround(trig)]++;
            }
            if (pass[2]) trigAfter[(int)std::lround(trig)]++;
            for (int s = 0; s < NS; ++s)
                if (pass[s]) { rows[s]++; if (b1) best[s]++; gp[s] = true; }
            if (pass[3]) rr[2 + 2 * k] += 1;
        }
        flush();
        t->ResetBranchAddresses();
        const char* sname[NS] = {"all", "isBest>=0", "vz>-500", "PT3"};
        for (int s = 0; s < NS; ++s)
            fprintf(o, "STEP %s %d %s %lld %lld %lld\n", tag[k], s, sname[s], rows[s], best[s], evts[s]);
        fprintf(o, "CHECK %s isBestNeg %lld groupsBestNot1 %lld vzNoVertex(<-999) %lld vzFail(<=-500) %lld\n",
                tag[k], isBestNeg, groupsBestNot1, vzNoVertex, vzLow);
        for (const auto& p : trigBefore) fprintf(o, "TRIG %s before %d %lld\n", tag[k], p.first, p.second);
        for (const auto& p : trigAfter) fprintf(o, "TRIG %s after %d %lld\n", tag[k], p.first, p.second);
        fh->cd();
        hVz->Write(); hVzB->Write();
    }
    fh->Close();
    delete fh;

    // ---- runs ----
    int nRunsEpEm = 0, nRunsUS = 0, nRunsUS3 = 0;
    FILE* fr = fopen(Form("%s/runs_%03d.txt", outDir.c_str(), day), "w");
    fprintf(fr, "# run EpEm_ID_rows us_all us_PT3 pp_all pp_PT3 mm_all mm_PT3 (day %03d)\n", day);
    for (const auto& p : runRows) {
        const auto& v = p.second;
        if (v[0] > 0) nRunsEpEm++;
        if (v[1] > 0) nRunsUS++;
        if (v[2] > 0) nRunsUS3++;
        fprintf(fr, "%d %lld %lld %lld %lld %lld %lld %lld\n", p.first, v[0], v[1], v[2], v[3], v[4], v[5], v[6]);
    }
    fclose(fr);
    fprintf(o, "RUNS %d EpEm_ID %d us_all %d us_PT3 %d\n", day, nRunsEpEm, nRunsUS, nRunsUS3);
    fprintf(o, "TIME %.1f\n", sw.RealTime());
    fclose(o);
    f->Close();
    delete f;
    printf("smp_pat_counts: day %03d done in %.1f s\n", day, sw.RealTime());
}
