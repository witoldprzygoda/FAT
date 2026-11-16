// fix2_x17count.C -- counts the S0 rows of the X17 counting window with the rounded edges (Part III, 0.5312-0.5566)
// and with the exact edges of Part II (simulation window x 0.5439/0.5460), M(ee) 15.5-24.5 MeV/c2. Read-only.
void fix2_x17count() {
  TFile* f = TFile::Open("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root");
  TTree* t = (TTree*)f->Get("us");
  const char* base = "isBest>=0 && vz>-500 && vz>-160 && cor_ok==1 && pip_isUsedVertex==1 && pim_isUsedVertex==1 && mee>=0.0155 && mee<0.0245";
  const double r = 0.5439 / 0.5460;
  printf("ratio %.8f exact edges %.6f - %.6f\n", r, 0.5333 * r, 0.5587 * r);
  Long64_t a = t->GetEntries(Form("%s && m4>=0.5312 && m4<0.5566", base));
  Long64_t b = t->GetEntries(Form("%s && m4>=%.8f && m4<%.8f", base, 0.5333 * r, 0.5587 * r));
  printf("S0 X17 window rows: rounded edges 0.5312-0.5566: %lld ; exact edges: %lld\n", a, b);
}
