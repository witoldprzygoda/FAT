// ecalacc_emcgeo_vae25.C -- package ecalacc (2026-11-15): ECAL geometry and cell status of the SIMULATION
// Reads EmcGeomPar and EmcCellStatusPar of the parameter file used by the gen4 DST production of the SMASH simulation
//   /lustre/hades/user/hadesdst/sub/feb22/sim/p4500p/gen4/smash/smash_lepton/dst/allParam_feb22_sim_run_17000_gen4_10032026.root
// (analysisDST.cc there: paramSource root, release FEB22SIM_v4a, refId 17000) and prints one line per module:
//   CELL sector cell status theta_c phi_c x_c y_c z_c  x1 y1 z1 ... x8 y8 z8   (lab, mm; corners of the module)
// Lab positions as in Hydra 2.7.6 emc/hemcclusterf.cc (module transform, then the sector lab transform).
// Unnamed macro, piped to ROOT on vae25 (nothing is written there):
//   (cat plots/all_report/ecalacc_emcgeo_vae25.C; echo ".q") | ssh -o RemoteCommand=none -o RequestTTY=no vae25 \
//     'source /lustre/hades/user/przygoda/pat3/profile.sh >/dev/null 2>&1; cd /tmp; root -l -b 2>&1' > emcgeo_vae25_stdout_v2.txt
//   grep "^CELL" emcgeo_vae25_stdout_v2.txt > /hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/emc_cells_sim_v2.txt
#include "emcdef.h"
{
  gSystem->Load("libHydra"); gSystem->Load("libEmc"); gSystem->Load("libDst");
  Hades* myHades = new Hades();
  Int_t mdcMods[6][4] = {{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1}};
  HDst::setupSpectrometer("feb22", mdcMods, "emc");
  TString par = "/lustre/hades/user/hadesdst/sub/feb22/sim/p4500p/gen4/smash/smash_lepton/dst/allParam_feb22_sim_run_17000_gen4_10032026.root";
  HDst::setupParameterSources("root", "", par, "FEB22SIM_v4a");
  HRuntimeDb* rtdb = gHades->getRuntimeDb();
  HEmcGeomPar* gp = (HEmcGeomPar*) rtdb->getContainer("EmcGeomPar");
  HEmcCellStatusPar* sp = (HEmcCellStatusPar*) rtdb->getContainer("EmcCellStatusPar");
  rtdb->initContainers(17000);
  printf("BEGINTABLE\n");
  for (Int_t s = 0; s < 6; s++) {
    HModGeomPar* m = gp->getModule(s);
    if (!m) { printf("NOSECTOR %d\n", s); continue; }
    HGeomTransform lab = m->getLabTransform();
    HGeomCompositeVolume* mod = m->getRefVolume();
    if (s == 0) { HGeomVector t = lab.getTransVector(); printf("SEC0TRANS %f %f %f\n", t.getX(), t.getY(), t.getZ()); }
    for (Int_t c = 0; c < emcMaxComponents; c++) {
      HGeomVolume* v = mod->getComponent(c);
      if (!v || v->getNumPoints() != 8) continue;
      HGeomVector p = v->getTransform().getTransVector();
      HGeomVector pl = lab.transFrom(p);
      Double_t th = TMath::ATan2(TMath::Sqrt(pl.getX()*pl.getX()+pl.getY()*pl.getY()), pl.getZ())*TMath::RadToDeg();
      Double_t ph = TMath::ATan2(pl.getY(), pl.getX())*TMath::RadToDeg(); if (ph < 0) ph += 360.;
      Int_t st = sp->getCellStatus(s, emc_cell2tile_map[c]);   // the status is indexed by the tile (position) number, as in hemcdigitizer.cc
      printf("CELL %d %d %d %.3f %.3f %.2f %.2f %.2f", s, c, st, th, ph, pl.getX(), pl.getY(), pl.getZ());
      for (Int_t i = 0; i < 8; i++) {
        HGeomVector q = *(v->getPoint(i));
        HGeomVector ql = lab.transFrom(v->getTransform().transFrom(q));
        printf(" %.2f %.2f %.2f", ql.getX(), ql.getY(), ql.getZ());
      }
      printf("\n");
    }
  }
  printf("ENDTABLE\n");
}
