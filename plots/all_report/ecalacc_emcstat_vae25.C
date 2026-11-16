// ecalacc_emcstat_vae25.C -- package ecalacc (2026-11-15): ECAL cell status of the DATA (feb22 gen4 real parameters)
// Reads EmcGeomPar and EmcCellStatusPar of /cvmfs/hadessoft.gsi.de/param/real/feb22/allParam_Feb22_gen4_07032026.root for
// eight runs (first runs of days 034, 039, 045, 047, 049, 053, 059, 067) and prints the modules per sector, the active ones
// and (first run) the list of status-0 modules. Unnamed macro, piped to ROOT on vae25 (nothing is written there):
//   (cat plots/all_report/ecalacc_emcstat_vae25.C; echo ".q") | ssh -o RemoteCommand=none -o RequestTTY=no vae25 \
//     'source /lustre/hades/user/przygoda/pat3/profile.sh >/dev/null 2>&1; cd /tmp; root -l -b 2>&1' > emcstat_real_vae25_stdout_v2.txt
#include "emcdef.h"
{
  gSystem->Load("libHydra"); gSystem->Load("libEmc"); gSystem->Load("libDst");
  Hades* myHades = new Hades();
  Int_t mdcMods[6][4] = {{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1}};
  HDst::setupSpectrometer("feb22", mdcMods, "emc");
  TString par = "/cvmfs/hadessoft.gsi.de/param/real/feb22/allParam_Feb22_gen4_07032026.root";
  HDst::setupParameterSources("root", "", par, "now");
  HRuntimeDb* rtdb = gHades->getRuntimeDb();
  HEmcGeomPar* gp = (HEmcGeomPar*) rtdb->getContainer("EmcGeomPar");
  HEmcCellStatusPar* sp = (HEmcCellStatusPar*) rtdb->getContainer("EmcCellStatusPar");
  Int_t runs[8] = {443850846, 444297309, 444793253, 444966037, 445138801, 445484442, 446061982, 446694022};
  for (Int_t ir = 0; ir < 8; ir++) {
    rtdb->initContainers(runs[ir]);
    Int_t nact[6] = {0,0,0,0,0,0}, ngeo[6] = {0,0,0,0,0,0};
    for (Int_t s = 0; s < 6; s++) {
      HModGeomPar* m = gp->getModule(s);
      if (!m) continue;
      HGeomCompositeVolume* mod = m->getRefVolume();
      for (Int_t c = 0; c < emcMaxComponents; c++) {
        HGeomVolume* v = mod->getComponent(c);
        if (!v || v->getNumPoints() != 8) continue;
        ngeo[s]++;
        if (sp->getCellStatus(s, emc_cell2tile_map[c]) == 1) nact[s]++;
      }
    }
    printf("RUN %d geo %d %d %d %d %d %d active %d %d %d %d %d %d\n", runs[ir], ngeo[0],ngeo[1],ngeo[2],ngeo[3],ngeo[4],ngeo[5], nact[0],nact[1],nact[2],nact[3],nact[4],nact[5]);
    if (ir == 0) {
      HModGeomPar* m = gp->getModule(1);
      if (m) { HGeomVector t = m->getLabTransform().getTransVector(); printf("SEC1TRANS %f %f %f\n", t.getX(), t.getY(), t.getZ()); }
      for (Int_t s = 0; s < 6; s++) for (Int_t c = 0; c < emcMaxComponents; c++) {
        HModGeomPar* mm = gp->getModule(s); if (!mm) continue;
        HGeomVolume* v = mm->getRefVolume()->getComponent(c);
        if (!v || v->getNumPoints() != 8) continue;
        if (sp->getCellStatus(s, emc_cell2tile_map[c]) != 1) printf("DEAD %d %d %d\n", runs[ir], s, c);
      }
    }
  }
}
