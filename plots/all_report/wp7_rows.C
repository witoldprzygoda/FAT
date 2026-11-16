// ============================================================================
// wp7_rows.C  --  WP7 (systematics, SPEC report_all): data rows for every selection variant
// ============================================================================
// Written by WP7 (prefix wp7), 2026-11-14.
// One pass over the data compact file (scan_v2). Every row (e+e-, e+e+, e-e-) with isBest >= 0, vz > -500, PT3 and
// vz > -170 mm is written to a small tree, in file order (event groups stay adjacent), with the columns that the
// selection variants of WP7 need:
//   k (0 e+e-, 1 e+e+, 2 e-e-), grp, kinematics m4, mpipi, mee, mm4, pt4, y4 (cor_ok == 0 rows: the *_fb values, as
//   wp0_fill.C ValsOf for the cor_ok fallback), vz, dz_pi, isBest, cor (cor_ok), vtx (both pions in the vertex fit),
//   per (alldefs period), fold, day.
// Selections are applied when the tree is read (wp7_lib.h, PassRow): S0, pre, dz40, PEAK, vz150, vz170, isBest == 1,
// cor_ok fallback, periods, folds. The check at the end compares the window counts with hists_v1.root.
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/wp7/wp7_rows.root (+ wp7_rows.txt)
// Run (from /home/przygoda/HADES/FAT/FAT; about 1 min, 3 GB):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp7", kTRUE)' 'plots/all_report/wp7_rows.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"
#include "hists_io.h"

void wp7_rows() {
    using namespace alldefs;
    const char* out = "/hdd3/przygoda/HADES/FAT_out/all_report/wp7";
    gSystem->mkdir(out, kTRUE);
    FILE* log = fopen(Form("%s/wp7_rows.txt", out), "w");
    scan::Data d;
    if (!scan::Load(d)) return;
    TFile* fo = TFile::Open(Form("%s/wp7_rows.root", out), "RECREATE");
    TTree* t = new TTree("rows", "WP7 rows: isBest >= 0, vz > -500, PT3, vz > -170; k 0 e+e-, 1 e+e+, 2 e-e-; file order; cor_ok == 0 rows carry the _fb kinematics");
    Char_t k, isBest, cor, vtx, per, fold; Short_t day; Int_t grp;
    Float_t m4, mpipi, mee, mm4, pt4, y4, vz, dz;
    t->Branch("k", &k, "k/B"); t->Branch("grp", &grp, "grp/I");
    t->Branch("m4", &m4, "m4/F"); t->Branch("mpipi", &mpipi, "mpipi/F"); t->Branch("mee", &mee, "mee/F");
    t->Branch("mm4", &mm4, "mm4/F"); t->Branch("pt4", &pt4, "pt4/F"); t->Branch("y4", &y4, "y4/F");
    t->Branch("vz", &vz, "vz/F"); t->Branch("dz", &dz, "dz/F");
    t->Branch("isBest", &isBest, "isBest/B"); t->Branch("cor", &cor, "cor/B"); t->Branch("vtx", &vtx, "vtx/B");
    t->Branch("per", &per, "per/B"); t->Branch("fold", &fold, "fold/B"); t->Branch("day", &day, "day/S");
    long n[3] = {0, 0, 0}, nS0[3] = {0, 0, 0};
    double nP = 0, nB = 0;
    for (int s = 0; s < 3; ++s)
        for (const auto& r : d.s[s]) {
            if (!PassBase(r, true) || !(r.vz > -170.f)) continue;
            const bool fb = !r.cor_ok;
            k = s; grp = r.grp;
            m4 = fb ? r.m4_fb : r.m4; mpipi = fb ? r.mpipi_fb : r.mpipi; mee = fb ? r.mee_fb : r.mee;
            mm4 = fb ? r.mm4_fb : r.mm4; pt4 = fb ? r.pt4_fb : r.pt4; y4 = fb ? r.y4_fb : r.y4;
            vz = r.vz; dz = r.dz_pi; isBest = r.isBest; cor = r.cor_ok ? 1 : 0; vtx = PiVtx(r) ? 1 : 0;
            per = PeriodOfDay(r.file_id); fold = r.fold & 1; day = r.file_id;
            t->Fill(); ++n[s];
            if (PassSel(kSelS0, r, r.cor_ok == 1, r.dz_pi, r.mm4)) {
                ++nS0[s];
                if (s == 0 && r.m4 >= kPlo && r.m4 < kPhi) nP += 1;
                if (s == 0 && r.m4 >= kBlo && r.m4 < kBhi) nB += 1;
            }
        }
    fo->cd(); t->Write(); fo->Close();
    TFile* fh = TFile::Open(kHistFile);
    TH1* nw = allio::Data(fh, "S0", "all", "N", "win");
    const char* line = Form("rows written: e+e- %ld, e+e+ %ld, e-e- %ld; S0 (cor_ok): %ld / %ld / %ld; S0 N in P %.0f (hists_v1 %.0f), in B %.0f (hists_v1 %.0f)\n",
                            n[0], n[1], n[2], nS0[0], nS0[1], nS0[2], nP, nw->GetBinContent(kWP + 1), nB, nw->GetBinContent(kWB + 1));
    printf("%s", line); fprintf(log, "%s", line);
    fclose(log);
}
