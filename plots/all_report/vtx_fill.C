// ==============================================================================================================
// vtx_fill.C -- copy (2026-11-14, integration of the event-vertex efficiency into report_all / report_exp) of
//   /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_verify/vv_fill.C
//   (study of the eta efficiency of the mandatory event-vertex cut vz > -500 mm; verifier ve_verify, own code).
// What it does: rows straight from the PAT files (data and simulation), WITHOUT any vertex requirement.
// Changes against the original: only the file and function names (vv_fill -> vtx_fill); the rows are written to the same
//   place as in the original run, /hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/rows/{d,s}_XX.root (re-running rewrites them).
// The original outputs are copied to docs/latex/tables_all/vtx_*.txt; the notes to docs/latex/notes_all/ve_*_note.md.
// Run (from /home/przygoda/HADES/FAT/FAT; vtx_fill.C first, then the analysis macros, at most 6 processes at once):
//   B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic", kTRUE)'
//   seq 0 28 | xargs -P 6 -I{} sh -c "root -l -b -q -e '$B' 'plots/all_report/vtx_fill.C+({},0)'"   # data, 29 day files
//   for i in 0 1; do root -l -b -q -e "$B" "plots/all_report/vtx_fill.C+($i,1)"; done                 # simulation (eta, std set)
// ==============================================================================================================
// vtx_fill.C -- ve_verify: independent row extraction from the PAT day files (data) or the simulated PAT files
// Usage: root -l -b -q 'vtx_fill.C+(ifile, isSim)'
//   data: ifile = 0..28 (line of lep_betanew.list); trees us/pp/mm; rows with isBest >= 0 and trigbit == 8192
//   sim : ifile = 0 (ppee_eta_all), 1 (ppee_std_all); tree us only; rows with isBest >= 0 and (trigbit & 8192) != 0
// Output: /hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/rows/{d,s}_XX.root, tree "t" (one row per stored combination,
// NO vertex requirement). Written from scratch for the verification (not derived from the ve_data / ve_sim_cause macros).
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TString.h"

namespace vv {
const double kMpi = 139.56995, kMe = 0.51099895, kD2R = M_PI / 180.0;
TLorentzVector P4(double p, double th, double ph, double m) {
  TLorentzVector v; v.SetXYZM(p * std::sin(th * kD2R) * std::cos(ph * kD2R), p * std::sin(th * kD2R) * std::sin(ph * kD2R),
                              p * std::cos(th * kD2R), m);
  return v;
}
// closest approach of two track lines; HYDRA line: base (-r sin phi, r cos phi, z), dir (theta, phi)
void Cross(double r1, double z1, double t1, double p1, double r2, double z2, double t2, double p2, float& zx, float& rx, float& dca) {
  t1 *= kD2R; p1 *= kD2R; t2 *= kD2R; p2 *= kD2R;
  TVector3 b1(-r1 * std::sin(p1), r1 * std::cos(p1), z1), b2(-r2 * std::sin(p2), r2 * std::cos(p2), z2);
  TVector3 d1(std::sin(t1) * std::cos(p1), std::sin(t1) * std::sin(p1), std::cos(t1));
  TVector3 d2(std::sin(t2) * std::cos(p2), std::sin(t2) * std::sin(p2), std::cos(t2));
  TVector3 w = b1 - b2; double b = d1.Dot(d2), d = d1.Dot(w), e = d2.Dot(w), den = 1 - b * b;
  if (den < 1e-12) { zx = -9999; rx = -1; dca = -1; return; }
  double s = (b * e - d) / den, t = (e - b * d) / den;
  TVector3 q1 = b1 + s * d1, q2 = b2 + t * d2, q = 0.5 * (q1 + q2);
  zx = q.Z(); rx = std::hypot(q.X(), q.Y()); dca = (q1 - q2).Mag();
}
}  // namespace vv

void vtx_fill(int ifile = 0, int isSim = 0) {
  using namespace vv;
  std::string path;
  if (!isSim) {
    std::ifstream in("/hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list");
    std::string l; int k = 0;
    while (std::getline(in, l)) { if (l.empty()) continue; if (k == ifile) path = l; ++k; }
  } else {
    path = ifile == 0 ? "/hdd3/przygoda/HADES/files/sim/ppee_eta_all.root" : "/hdd3/przygoda/HADES/files/sim/ppee_std_all.root";
  }
  if (path.empty()) { printf("no file %d\n", ifile); return; }
  int day = -1;
  if (!isSim) { size_t p = path.find("lep_"); day = std::atoi(path.substr(p + 4, 3).c_str()); }
  TFile* f = TFile::Open(path.c_str());
  TString out = Form("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/rows/%s_%02d.root", isSim ? "s" : "d", ifile);
  TFile* fo = new TFile(out, "RECREATE");
  TTree* to = new TTree("t", "ve_verify rows");
  // output columns
  Int_t o_smp, o_day, o_run, o_evt, o_grp, o_best, o_start, o_vcand, o_cor, o_ng, o_trig;
  Float_t o_vz, o_vx, o_vy, o_vcl, o_m4, o_mee, o_mpp, o_oa, o_meeg, o_m5c, o_zx, o_rx, o_dca, o_zxe, o_rxe;
  Float_t o_z[4], o_r[4], o_th[4];
  // sim truth: per track id, parent, grandparent, process, geninfo, geninfo1, geninfo2, vertex
  Int_t o_tid[4], o_tpar[4], o_tgp[4], o_tproc[4];
  Float_t o_g0[4], o_g1[4], o_g2[4], o_tvz[4], o_tvx[4], o_tvy[4];
  to->Branch("smp", &o_smp); to->Branch("day", &o_day); to->Branch("run", &o_run); to->Branch("evt", &o_evt);
  to->Branch("grp", &o_grp); to->Branch("best", &o_best); to->Branch("start", &o_start); to->Branch("vcand", &o_vcand);
  to->Branch("cor", &o_cor); to->Branch("ng", &o_ng); to->Branch("trig", &o_trig);
  to->Branch("vz", &o_vz); to->Branch("vx", &o_vx); to->Branch("vy", &o_vy); to->Branch("vcl", &o_vcl);
  to->Branch("m4", &o_m4); to->Branch("mee", &o_mee); to->Branch("mpp", &o_mpp); to->Branch("oa", &o_oa);
  to->Branch("meeg", &o_meeg); to->Branch("m5c", &o_m5c);
  to->Branch("zx", &o_zx); to->Branch("rx", &o_rx); to->Branch("dca", &o_dca); to->Branch("zxe", &o_zxe); to->Branch("rxe", &o_rxe);
  to->Branch("z", o_z, "z[4]/F"); to->Branch("r", o_r, "r[4]/F"); to->Branch("th", o_th, "th[4]/F");
  if (isSim) {
    to->Branch("tid", o_tid, "tid[4]/I"); to->Branch("tpar", o_tpar, "tpar[4]/I"); to->Branch("tgp", o_tgp, "tgp[4]/I");
    to->Branch("tproc", o_tproc, "tproc[4]/I"); to->Branch("g0", o_g0, "g0[4]/F"); to->Branch("g1", o_g1, "g1[4]/F");
    to->Branch("g2", o_g2, "g2[4]/F"); to->Branch("tvz", o_tvz, "tvz[4]/F"); to->Branch("tvx", o_tvx, "tvx[4]/F");
    to->Branch("tvy", o_tvy, "tvy[4]/F");
  }
  const char* trees[3] = {"PipPimEpEm_ID", "PipPimEpEp_ID", "PipPimEmEm_ID"};
  const char* l1[3] = {"ep", "ep1", "em1"};
  const char* l2[3] = {"em", "ep2", "em2"};
  long long nIn[3] = {0, 0, 0}, nOut[3] = {0, 0, 0};
  int grp = 0;
  for (int s = 0; s < (isSim ? 1 : 3); ++s) {
    TTree* t = (TTree*)f->Get(trees[s]);
    if (!t) continue;
    t->SetBranchStatus("*", 0);
    std::vector<std::string> names;
    auto B = [&](const std::string& n, Float_t* a) { t->SetBranchStatus(n.c_str(), 1); t->SetBranchAddress(n.c_str(), a); };
    Float_t run, evt, best, trig, start, vcand, vx, vy, vz, vcl, nmult;
    B("runnumber", &run); B("event", &evt); B("isBest", &best); B("trigbit", &trig); B("evtStart", &start);
    B("evtVertCand", &vcand); B("eVertReco_x", &vx); B("eVertReco_y", &vy); B("eVertReco_z", &vz); B("eVertClust_z", &vcl);
    B("neutr_mult", &nmult);
    const std::string pre[4] = {l1[s], l2[s], "pip", "pim"};
    const std::string corr[4] = {"_p_corr_ep", "_p_corr_em", "_p_corr_pip", "_p_corr_pim"};
    Float_t pc[4], th[4], ph[4], zz[4], rr[4];
    Float_t sid[4], spar[4], sgp[4], sproc[4], g0[4], g1[4], g2[4], tvx[4], tvy[4], tvz[4];
    for (int k = 0; k < 4; ++k) {
      B(pre[k] + corr[k], &pc[k]); B(pre[k] + "_theta", &th[k]); B(pre[k] + "_phi", &ph[k]);
      B(pre[k] + "_z", &zz[k]); B(pre[k] + "_r", &rr[k]);
      if (isSim) {
        B(pre[k] + "_sim_id", &sid[k]); B(pre[k] + "_sim_parentid", &spar[k]); B(pre[k] + "_sim_grandparentid", &sgp[k]);
        B(pre[k] + "_sim_processid", &sproc[k]); B(pre[k] + "_sim_geninfo", &g0[k]); B(pre[k] + "_sim_geninfo1", &g1[k]);
        B(pre[k] + "_sim_geninfo2", &g2[k]); B(pre[k] + "_sim_vertexx", &tvx[k]); B(pre[k] + "_sim_vertexy", &tvy[k]);
        B(pre[k] + "_sim_vertexz", &tvz[k]);
      }
    }
    Float_t npid[5], nbeta[5], ne[5], nth[5], nph[5];
    for (int k = 0; k < 5; ++k) {
      B(Form("neutr_pid_%d", k + 1), &npid[k]); B(Form("neutr_beta_%d", k + 1), &nbeta[k]);
      B(Form("neutr_cluster_energy_%d", k + 1), &ne[k]); B(Form("neutr_cluster_theta_%d", k + 1), &nth[k]);
      B(Form("neutr_cluster_phi_%d", k + 1), &nph[k]);
    }
    const Long64_t n = t->GetEntries();
    float lrun = -1, levt = -1;
    for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      // event group: adjacent rows with equal (run, event), counted on ALL rows of the tree
      if (run != lrun || evt != levt) { ++grp; lrun = run; levt = evt; }
      if (!(best >= 0)) continue;
      const int tb = (int)std::lround(trig);
      if (isSim ? ((tb & 8192) == 0) : (tb != 8192)) continue;
      ++nIn[s];
      o_smp = s; o_day = day; o_run = (int)std::lround(run); o_evt = (int)std::lround(evt); o_grp = grp;
      o_best = (int)std::lround(best); o_start = (int)std::lround(start); o_vcand = (int)std::lround(vcand); o_trig = tb;
      o_vz = vz; o_vx = vx; o_vy = vy; o_vcl = vcl;
      bool ok = true;
      for (int k = 0; k < 4; ++k) if (!std::isfinite(pc[k])) ok = false;
      o_cor = ok;
      TLorentzVector e1 = P4(pc[0], th[0], ph[0], kMe), e2 = P4(pc[1], th[1], ph[1], kMe);
      TLorentzVector pp = P4(pc[2], th[2], ph[2], kMpi), pm = P4(pc[3], th[3], ph[3], kMpi);
      TLorentzVector ee = e1 + e2, pipi = pp + pm, v4 = ee + pipi;
      o_m4 = ok ? v4.M() / 1000. : NAN; o_mee = ok ? ee.M() / 1000. : NAN; o_mpp = ok ? pipi.M() / 1000. : NAN;
      o_oa = e1.Vect().Angle(e2.Vect()) / kD2R;
      // photon: neutr_pid == 1, 0.8 <= beta <= 1.2, cluster energy > 100 MeV; 3 most energetic; M(ee g) closest to 0.140
      o_ng = 0; o_meeg = -1; o_m5c = -1;
      {
        int nm = std::min(5, std::max(0, (int)std::lround(nmult)));
        std::vector<std::pair<float, int>> g;
        for (int k = 0; k < nm; ++k)
          if (std::lround(npid[k]) == 1 && nbeta[k] >= 0.8f && nbeta[k] <= 1.2f && ne[k] > 100.f) g.push_back({ne[k], k});
        std::sort(g.begin(), g.end(), [](const std::pair<float, int>& a, const std::pair<float, int>& b) { return a.first > b.first; });
        o_ng = g.size();
        double bestd = 1e9;
        for (size_t j = 0; j < g.size() && j < 3 && ok; ++j) {
          const int k = g[j].second;
          TLorentzVector vg = P4(ne[k], nth[k], nph[k], 0.);
          const double m = (ee + vg).M() / 1000.;
          if (std::fabs(m - 0.140) < bestd) {
            bestd = std::fabs(m - 0.140); o_meeg = m;
            const double mee2 = ee.M2() / 1e6, mg2 = m * m, mpi0 = 0.1349768;
            const double fsc = (mg2 - mee2) > 0 ? (mpi0 * mpi0 - mee2) / (mg2 - mee2) : 1.;
            o_m5c = (v4 + fsc * vg).M() / 1000.;
          }
        }
      }
      for (int k = 0; k < 4; ++k) { o_z[k] = zz[k]; o_r[k] = rr[k]; o_th[k] = th[k]; }
      Cross(rr[2], zz[2], th[2], ph[2], rr[3], zz[3], th[3], ph[3], o_zx, o_rx, o_dca);
      float dd; Cross(rr[0], zz[0], th[0], ph[0], rr[1], zz[1], th[1], ph[1], o_zxe, o_rxe, dd);
      if (isSim)
        for (int k = 0; k < 4; ++k) {
          o_tid[k] = std::lround(sid[k]); o_tpar[k] = std::lround(spar[k]); o_tgp[k] = std::lround(sgp[k]);
          o_tproc[k] = std::lround(sproc[k]); o_g0[k] = g0[k]; o_g1[k] = g1[k]; o_g2[k] = g2[k];
          o_tvx[k] = tvx[k]; o_tvy[k] = tvy[k]; o_tvz[k] = tvz[k];
        }
      to->Fill(); ++nOut[s];
    }
  }
  fo->cd(); to->Write(); fo->Close();
  printf("file %d day %d %s : rows us %lld pp %lld mm %lld\n", ifile, day, path.c_str(), nIn[0], nIn[1], nIn[2]);
}
