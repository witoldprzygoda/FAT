// ==============================================================================================================
// vtx_load.h -- copy (2026-11-14, integration of the event-vertex efficiency into report_all / report_exp) of
//   /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_verify/vv_load.h
//   (study of the eta efficiency of the mandatory event-vertex cut vz > -500 mm; verifier ve_verify, own code).
// What it does: row loader.
// Changes against the original: only the file name and the include names (vv_* -> vtx_*); the code is unchanged.
// The original outputs are copied to docs/latex/tables_all/vtx_*.txt; the notes to docs/latex/notes_all/ve_*_note.md.
// Run (from /home/przygoda/HADES/FAT/FAT; vtx_fill.C first, then the analysis macros, at most 6 processes at once):
//   B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic", kTRUE)'
//   (header; reads the rows written by vtx_fill.C)
// ==============================================================================================================
std::vector<R> Load(const char* pat, bool sim = false) {
  TChain ch("t"); ch.Add(pat);
  Int_t smp, day, grp, best, start, vcand, cor, ng; Float_t vz, m4, mee, mpp, oa, meeg, m5c, zx, rx, dca, zxe, rxe, z[4], r[4];
  ch.SetBranchStatus("*", 0);
  auto BI = [&](const char* n, Int_t* a) { ch.SetBranchStatus(n, 1); ch.SetBranchAddress(n, a); };
  auto BF = [&](const char* n, Float_t* a) { ch.SetBranchStatus(n, 1); ch.SetBranchAddress(n, a); };
  BI("smp", &smp); BI("day", &day); BI("grp", &grp); BI("best", &best); BI("start", &start); BI("vcand", &vcand);
  BI("cor", &cor); BI("ng", &ng);
  BF("vz", &vz); BF("m4", &m4); BF("mee", &mee); BF("mpp", &mpp); BF("oa", &oa); BF("meeg", &meeg); BF("m5c", &m5c);
  BF("zx", &zx); BF("rx", &rx); BF("dca", &dca); BF("zxe", &zxe); BF("rxe", &rxe); BF("z", z); BF("r", r);
  std::vector<R> v; v.reserve(ch.GetEntries());
  for (Long64_t i = 0; i < ch.GetEntries(); ++i) {
    ch.GetEntry(i);
    R q; q.file = ch.GetTreeNumber(); q.grp = grp; q.smp = smp; q.day = day; q.best = best; q.start = start; q.vcand = vcand;
    q.ng = ng; q.vz = vz; q.m4 = m4; q.mee = mee; q.mpp = mpp; q.oa = oa; q.meeg = meeg; q.m5c = m5c; q.zx = zx; q.rx = rx;
    q.dca = dca; q.zxe = zxe; q.rxe = rxe; q.cor = cor;
    float zs[4] = {z[0], z[1], z[2], z[3]}; std::sort(zs, zs + 4); q.zp = 0.5f * (zs[1] + zs[2]);
    int nw = 0; for (int k = 0; k < 4; ++k) if (z[k] > -170 && z[k] < -30 && std::fabs(r[k]) < 20) ++nw; q.nwin = nw;
    v.push_back(q);
  }
  return v;
}
