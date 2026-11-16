// vecalacc_emcstat_vae25.C -- INDEPENDENT VERIFIER of package ecalacc (2026-11-15): raw ECAL status array of the DATA
// Prints the raw 6 x 163 EmcCellStatusPar array of the real feb22 gen4 parameters for runs spread over the beam time
// (first runs of days 034, 045, 053, 067: periods I, I, IIB, III).
// Run:
//   (cat plots/all_report/vecalacc_emcstat_vae25.C; echo ".q") | ssh -o RemoteCommand=none -o RequestTTY=no vae25 \
//     'source /lustre/hades/user/przygoda/pat3/profile.sh >/dev/null 2>&1; cd /tmp; root -l -b 2>&1' \
//     > /hdd3/przygoda/HADES/FAT_out/all_report/vecalacc/vemcstat_stdout.txt
{
  gSystem->Load("libHydra"); gSystem->Load("libEmc"); gSystem->Load("libDst");
  Hades* myHades = new Hades();
  Int_t mdcMods[6][4] = {{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1}};
  HDst::setupSpectrometer("feb22", mdcMods, "emc");
  HDst::setupParameterSources("root", "", "/cvmfs/hadessoft.gsi.de/param/real/feb22/allParam_Feb22_gen4_07032026.root", "now");
  HRuntimeDb* rtdb = gHades->getRuntimeDb();
  HEmcCellStatusPar* sp = (HEmcCellStatusPar*) rtdb->getContainer("EmcCellStatusPar");
  Int_t runs[4] = {443850846, 444793253, 445484442, 446694022};
  for (Int_t ir = 0; ir < 4; ir++) {
    rtdb->initContainers(runs[ir]);
    printf("RAWSTATUS_REAL %d", runs[ir]);
    for (Int_t i = 0; i < 6 * 163; i++) printf(" %d", sp->getCellStatus(i / 163, i % 163));
    printf("\n");
  }
}
