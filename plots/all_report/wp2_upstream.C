// ============================================================================
// wp2_upstream.C  --  WP2: vertex-z spectra of PAT rows in full-target and empty-target runs (flux normaliser)
// ============================================================================
// Plain words: reactions in the beam-line material upstream of the target happen at the same rate per beam
// particle whether the target cell is full or empty. Rows with a vertex far upstream therefore measure the
// recorded beam flux in both kinds of runs (with the same trigger and the same dead time as the yields).
// This macro fills the vertex-z spectra (5 mm bins, -700..+100 mm) of PT3 rows (trigbit 8192) of
//   EpEm_ID (isBest == 1)          -- high-statistics normaliser
//   PipPimEpEm_ID (isBest >= 0), PipPimEpEp_ID, PipPimEmEm_ID   -- the channel itself, before FAT
// for the full-target day files /hdd3/przygoda/HADES/files/exp/gen4/lep_<day>.root and the empty-target PAT outputs
// /hdd3/przygoda/HADES/FAT_out/all_report/wp2/empty/pat/lep_empty_<day>.root, days 034, 047, 065.
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_upstream.root (histograms h_<tree>_<full|empty>_<day>)
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp2", kTRUE)' \
//           'plots/all_report/wp2_upstream.C+O'
// ============================================================================
#include <cstdio>
#include "TFile.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TTree.h"

void wp2_upstream() {
    TStopwatch sw;
    const int days[3] = {34, 47, 65};
    const char* trees[4] = {"EpEm_ID", "PipPimEpEm_ID", "PipPimEpEp_ID", "PipPimEmEm_ID"};
    TFile* fo = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_upstream.root", "RECREATE");
    for (int s = 0; s < 2; ++s) {
        for (int d : days) {
            const TString fn = s == 0 ? Form("/hdd3/przygoda/HADES/files/exp/gen4/lep_%03d.root", d)
                                      : Form("/hdd3/przygoda/HADES/FAT_out/all_report/wp2/empty/pat/lep_empty_%03d.root", d);
            TFile* f = TFile::Open(fn, "READ");
            if (!f || f->IsZombie()) { printf("wp2_upstream: cannot open %s\n", fn.Data()); continue; }
            for (int k = 0; k < 4; ++k) {
                TTree* t = (TTree*)f->Get(trees[k]);
                if (!t) continue;
                Float_t vz, tb, isb;
                t->SetBranchStatus("*", 0);
                t->SetBranchStatus("eVertReco_z", 1); t->SetBranchStatus("trigbit", 1); t->SetBranchStatus("isBest", 1);
                t->SetBranchAddress("eVertReco_z", &vz); t->SetBranchAddress("trigbit", &tb); t->SetBranchAddress("isBest", &isb);
                fo->cd();
                TH1D* h = new TH1D(Form("h_%s_%s_%03d", trees[k], s == 0 ? "full" : "empty", d),
                                   Form("%s, %s target, day %03d, PT3, %s;v_{z} [mm];rows / 5 mm", trees[k], s == 0 ? "full" : "empty", d,
                                        k == 0 ? "isBest == 1" : "isBest >= 0"), 160, -700, 100);
                const Long64_t n = t->GetEntries();
                for (Long64_t i = 0; i < n; ++i) {
                    t->GetEntry(i);
                    if ((int)tb != 8192) continue;
                    if (k == 0 ? isb != 1.f : isb < 0.f) continue;
                    h->Fill(vz);
                }
                h->Write();
                printf("wp2_upstream: %-14s %s %03d: %lld rows, PT3 filled %.0f (%.1f s)\n", trees[k], s == 0 ? "full " : "empty", d, n,
                       h->GetEntries(), sw.RealTime());
                sw.Continue();
            }
            f->Close();
            delete f;
        }
    }
    fo->Close();
}
