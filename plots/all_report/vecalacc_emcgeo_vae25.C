// vecalacc_emcgeo_vae25.C -- INDEPENDENT VERIFIER of package ecalacc (2026-11-15): own dump of the ECAL geometry and status
// Reads EmcGeomPar and EmcCellStatusPar of the gen4 SMASH simulation parameter file (identical md5 b67f9133... in
// smash_standard/dst and smash_lepton/dst) and, for comparison, of the real feb22 gen4 parameter file at one run.
// Differences to ecalacc_emcgeo_vae25.C (on purpose):
//   - the status position is taken from HEmcDetector::getPositionFromCell (the call used by hemcdigitizer.cc), not from
//     the emc_cell2tile_map table; the raw 6 x 163 status array is printed as well;
//   - the front face is identified from the LOCAL z of the 8 box points (points with the smaller local z), not by the
//     distance to the origin; the local points are printed for one module.
// Output (stdout only, nothing is written on vae25):
//   VCELL sector cell pos status_sim  xc yc zc  (front face: 4 x (x y z))  (back face: 4 x (x y z))
//   (the real-data status is read by vecalacc_emcstat_vae25.C)
// Run:
//   (cat plots/all_report/vecalacc_emcgeo_vae25.C; echo ".q") | ssh -o RemoteCommand=none -o RequestTTY=no vae25 \
//     'source /lustre/hades/user/przygoda/pat3/profile.sh >/dev/null 2>&1; cd /tmp; root -l -b 2>&1' \
//     > /hdd3/przygoda/HADES/FAT_out/all_report/vecalacc/vemcgeo_stdout.txt
{
  gSystem->Load("libHydra"); gSystem->Load("libEmc"); gSystem->Load("libDst");
  Hades* myHades = new Hades();
  Int_t mdcMods[6][4] = {{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1}};
  HDst::setupSpectrometer("feb22", mdcMods, "emc");
  TString par = "/lustre/hades/user/hadesdst/sub/feb22/sim/p4500p/gen4/smash/smash_standard/dst/allParam_feb22_sim_run_17000_gen4_10032026.root";
  HDst::setupParameterSources("root", "", par, "FEB22SIM_v4a");
  HRuntimeDb* rtdb = gHades->getRuntimeDb();
  HEmcGeomPar* gp = (HEmcGeomPar*) rtdb->getContainer("EmcGeomPar");
  HEmcCellStatusPar* sp = (HEmcCellStatusPar*) rtdb->getContainer("EmcCellStatusPar");
  rtdb->initContainers(17000);
  Int_t stSim[6][255];
  for (Int_t s = 0; s < 6; s++) for (Int_t c = 0; c < 255; c++) {
    Int_t pos = HEmcDetector::getPositionFromCell(c);
    stSim[s][c] = pos >= 0 ? sp->getCellStatus(s, pos) : -9;
  }
  printf("RAWSTATUS_SIM");
  for (Int_t i = 0; i < 6 * 163; i++) printf(" %d", sp->getCellStatus(i / 163, i % 163));
  printf("\n");
  // print the local points of one module to see the box convention
  {
    HGeomVolume* v = gp->getModule(2)->getRefVolume()->getComponent(HEmcDetector::getPositionFromCell(0) >= 0 ? 6 : 6);
    if (v) for (Int_t i = 0; i < v->getNumPoints(); i++) { HGeomVector q = *(v->getPoint(i)); printf("LOCALPT %d %.3f %.3f %.3f\n", i, q.getX(), q.getY(), q.getZ()); }
    if (v) { HGeomVector t = v->getTransform().getTransVector(); printf("LOCALTRANS %.3f %.3f %.3f\n", t.getX(), t.getY(), t.getZ()); }
  }
  // geometry in the lab
  struct C { Int_t s, c, pos, st; Double_t x[8][3], ctr[3], lz[8]; };
  std::vector<C> cells;
  for (Int_t s = 0; s < 6; s++) {
    HModGeomPar* m = gp->getModule(s);
    if (!m) { printf("NOSECTOR %d\n", s); continue; }
    HGeomTransform lab = m->getLabTransform();
    HGeomVector t0 = lab.getTransVector();
    printf("SECTRANS %d %.3f %.3f %.3f\n", s, t0.getX(), t0.getY(), t0.getZ());
    HGeomCompositeVolume* mod = m->getRefVolume();
    for (Int_t c = 0; c < 255; c++) {
      HGeomVolume* v = mod->getComponent(c);
      if (!v || v->getNumPoints() != 8) continue;
      C q; q.s = s; q.c = c; q.pos = HEmcDetector::getPositionFromCell(c); q.st = stSim[s][c];
      HGeomVector p = lab.transFrom(v->getTransform().getTransVector());
      q.ctr[0] = p.getX(); q.ctr[1] = p.getY(); q.ctr[2] = p.getZ();
      for (Int_t i = 0; i < 8; i++) {
        HGeomVector loc = *(v->getPoint(i));
        q.lz[i] = loc.getZ();
        HGeomVector g = lab.transFrom(v->getTransform().transFrom(loc));
        q.x[i][0] = g.getX(); q.x[i][1] = g.getY(); q.x[i][2] = g.getZ();
      }
      cells.push_back(q);
    }
  }
  for (auto& q : cells) {
    // front face = the four points with the smaller LOCAL z (box: points 0-3 at -dz, 4-7 at +dz in Hydra's BOX convention)
    std::vector<Int_t> o = {0,1,2,3,4,5,6,7};
    std::sort(o.begin(), o.end(), [&](Int_t a, Int_t b) { return q.lz[a] < q.lz[b]; });
    printf("VCELL %d %d %d %d %.3f %.3f %.3f", q.s, q.c, q.pos, q.st, q.ctr[0], q.ctr[1], q.ctr[2]);
    for (Int_t k = 0; k < 8; k++) printf(" %.3f %.3f %.3f", q.x[o[k]][0], q.x[o[k]][1], q.x[o[k]][2]);
    printf("\n");
  }
  printf("NCELLS %zu\n", cells.size());
}
