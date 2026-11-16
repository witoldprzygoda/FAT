// ==============================================================================================================
// vtx_elsum.C -- copy (2026-11-14, integration of the event-vertex efficiency into report_all / report_exp) of
//   /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_sim_cause/elsum.C
//   (study of the eta efficiency of the mandatory event-vertex cut vz > -500 mm; package ve_sim_cause).
// What it does: sums the elastic excess: target cell (vertex) against upstream (no vertex); result 0.028 +- 0.003.
// Changes against the original: only the file and function names (-> vtx_*); the code is unchanged (it writes the ROOT file
//   given as argument; the summed result of the original run is ve_sim_cause/elsum_all.txt, quoted in notes_all/ve_sim_cause_note.md).
// The original outputs are copied to docs/latex/tables_all/vtx_*.txt; the notes to docs/latex/notes_all/ve_*_note.md.
// Run (from /home/przygoda/HADES/FAT/FAT; vtx_fill.C first, then the analysis macros, at most 6 processes at once):
//   B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic", kTRUE)'
//   root -l -b -q "plots/all_report/vtx_elsum.C(\"O/elastic_full_036_p1.root,O/elastic_full_049_p1.root,O/elastic_full_062_p1.root\")"   # O = /hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_sim_cause
// ==============================================================================================================
// vtx_elsum.C -- elastic content per 2-track event: target (vertex, zbar in [-160,-90)) vs upstream (no vertex, zbar < -170)
void vtx_elsum(const char* files) {
  TObjArray* fl = TString(files).Tokenize(",");
  double nT = 0, nU = 0, eT = 0, eU = 0, sT = 0, sU = 0;
  for (int i = 0; i < fl->GetEntries(); ++i) {
    TFile* f = TFile::Open(((TObjString*)fl->At(i))->GetString());
    TH1* v = (TH1*)f->Get("zbar_vtx"); TH1* n = (TH1*)f->Get("zbar_nov");
    nT += v->Integral(v->FindBin(-159.9), v->FindBin(-90.1)); nU += n->Integral(n->FindBin(-699.9), n->FindBin(-170.1));
    auto ex = [&](const char* nm, double& e, double& s) { TH1* h = (TH1*)f->Get(nm);
      double pk = h->Integral(h->FindBin(0.2543 + 1e-6), h->FindBin(0.3343 - 1e-6));
      double sb = 0.5 * (h->Integral(h->FindBin(0.15 + 1e-6), h->FindBin(0.23 - 1e-6)) + h->Integral(h->FindBin(0.36 + 1e-6), h->FindBin(0.44 - 1e-6)));
      e += pk - sb; s += pk + 0.5 * sb; };
    ex("tt_vtx_z5", eT, sT);
    for (int z = 0; z < 4; ++z) ex(Form("tt_nov_z%d", z), eU, sU);
  }
  printf("2-track events: target (vertex, zbar -160..-90) %.0f, upstream (no vertex, zbar < -170) %.0f\n", nT, nU);
  printf("elastic excess: target %.0f +- %.0f (%.4f per event), upstream %.0f +- %.0f (%.4f per event); ratio upstream/target %.3f +- %.3f\n",
         eT, sqrt(sT), eT / nT, eU, sqrt(sU), eU / nU, (eU / nU) / (eT / nT), sqrt(sU) / nU / (eT / nT));
}
