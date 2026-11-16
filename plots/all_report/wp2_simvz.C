// ============================================================================
// wp2_simvz.C  --  WP2: vertex-z shape of the simulated LH2 rows (template for the vz fit of the window share)
// ============================================================================
// SMASH std + lep + eta, owned rows (mode A), w x V3; selection = S0 WITHOUT the vz > -160 cut (isBest >= 0,
// vz > -500, PT3, cor_ok, both pions in the vertex fit). The simulation contains only the LH2 cell (no window
// material). Histograms (2 mm, -500..100) of the reconstructed vz for T12 (all M4, and the windows B and P),
// and of the true vertex of the pi+ (pip_sim_vz) for T12.
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_simvz.root
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp2", kTRUE)' 'plots/all_report/wp2_simvz.C+O'
// ============================================================================
#include <cstdio>
#include "TFile.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "defs.h"
#include "sim_eta_physics.h"

void wp2_simvz() {
    using namespace alldefs;
    TStopwatch sw;
    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    TFile* fo = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_simvz.root", "RECREATE");
    TH1D* h[4];
    const char* nm[4] = {"simvz_T12", "simvz_T12_B", "simvz_T12_P", "simvz_true_T12"};
    for (int k = 0; k < 4; ++k) { h[k] = new TH1D(nm[k], ";v_{z} [mm];#sigma [nb] / 2 mm", 300, -500, 100); h[k]->Sumw2(); }
    const char* files[3] = {simeng::kCompact, simeng::kCompact, simeng::kCompactEta};
    const char* trees[3] = {"std_us", "lep_us", "eta_us"};
    for (int t = 0; t < 3; ++t) {
        TFile* f = TFile::Open(files[t], "READ");
        TTree* tr = (TTree*)f->Get(trees[t]);
        simeng::Row r;
        simeng::BindRow(tr, r, false);
        const Long64_t n = tr->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            tr->GetEntry(i);
            if (r.ownA != 1) continue;
            if (!PassSel(kSelNone, r, r.cor_ok == 1, 0.f, r.mm4) || !PiVtx(r)) continue;
            int cn[6];
            const int ncn = CompsOf(r, false, cn);
            bool t12 = false;
            for (int q = 0; q < ncn; ++q) if (cn[q] == kT12) t12 = true;
            if (!t12) continue;
            const double w = r.w * simphys::RowFactor(r, kPhysNom) / simeng::kLsim;
            const double m4p = M4sim(r.m4, r.m4_true);
            h[0]->Fill(r.vz, w);
            if (m4p >= kBlo && m4p < kBhi) h[1]->Fill(r.vz, w);
            if (m4p >= kPlo && m4p < kPhi) h[2]->Fill(r.vz, w);
            h[3]->Fill(r.pip_sim_vz, w);
        }
        tr->ResetBranchAddresses();
        f->Close();
        printf("wp2_simvz: %s done (%.1f s)\n", trees[t], sw.RealTime());
        sw.Continue();
    }
    fo->cd();
    for (int k = 0; k < 4; ++k) h[k]->Write();
    const double all = h[0]->Integral(), below = h[0]->Integral(1, h[0]->FindBin(-160.001));
    printf("wp2_simvz: T12 %.3f nb, below vz = -160: %.4f %%, in [-200,-165): %.4f %%\n", all, 100 * below / all,
           100 * h[0]->Integral(h[0]->FindBin(-199.999), h[0]->FindBin(-165.001)) / all);
    fo->Close();
}
