// ============================================================================
// verify_spec_vz500.C  --  independent check of the PAT-level cut flow and of
//                          the pass fraction of the FAT cut "event vertex
//                          z > -500 mm" (chapter "spectra and kinematics")
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Written by the verifier of the spec chapter; does not use spec_common.h or
// ScanEngine.h. For the files with index i % nparts == part of
// /hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list it reads the trees
// PipPimEpEm_ID, PipPimEpEp_ID, PipPimEmEm_ID and counts rows (isBest >= 0) at
//   PAT -> eVertReco_z > -500 -> + trigbit == 8192 -> + cor_ok -> + vz > -160,
// the rows failing z > -500 that sit exactly at -1000 and their evtVertCand,
// and the z > -500 pass fraction of cor_ok rows in the data peak window
// 0.53364-0.55416 and the data bump window 0.4550-0.5225 GeV/c2, with M4 built
// with TLorentzVector from the corrected momenta as FAT reads them (lepton 1:
// <p1>_p_corr_ep, lepton 2: <p2>_p_corr_em, pions <pi>_p_corr_<pi>).
// Output: one line per sample "VZ500 k n0 n1 n2 n3 n4 nfail nm1000 nfailCand0
//         tPk pPk tBp pBp" on stdout (summed by the verifier with awk).
// Run (from /home/przygoda/HADES/FAT/FAT), e.g. six parts:
//   for k in 0 1 2 3 4 5; do root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_spec_verify", kTRUE)' "plots/exp_report/verify_spec_vz500.C+($k,6)" > /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/spec_verify/vz500_$k.log 2>&1 & done; wait
//   grep -h ^VZ500 /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/spec_verify/vz500_*.log | awk '{for(i=3;i<=NF;i++) s[$2,i]+=$i} END{for(k=0;k<3;k++){printf "%d",k; for(i=3;i<=14;i++) printf " %d", s[k,i]; printf "\n"}}'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"

void verify_spec_vz500(int part = 0, int nparts = 1) {
    std::vector<std::string> files;
    {
        std::ifstream in("/hdd3/przygoda/HADES/files/exp/gen4/lep_betanew.list");
        std::string s;
        while (std::getline(in, s)) if (!s.empty() && s[0] != '#') files.push_back(s);
    }
    const char* tree[3] = {"PipPimEpEm_ID", "PipPimEpEp_ID", "PipPimEmEm_ID"};
    const char* p1[3] = {"ep", "ep1", "em1"};
    const char* p2[3] = {"em", "ep2", "em2"};
    const double me = 0.51099895, mpi = 139.57039;
    const double pkLo = 0.5439 - 1.9 * 0.0054, pkHi = 0.5439 + 1.9 * 0.0054, bpLo = 0.455, bpHi = 0.5225;
    long long c[3][12] = {{0}};
    for (size_t fi = 0; fi < files.size(); ++fi) {
        if ((int)(fi % nparts) != part) continue;
        TFile* f = TFile::Open(files[fi].c_str());
        for (int k = 0; k < 3; ++k) {
            TTree* t = (TTree*)f->Get(tree[k]);
            t->SetBranchStatus("*", 0);
            std::vector<std::string> nm = {"isBest", "trigbit", "eVertReco_z", "evtVertCand",
                std::string(p1[k]) + "_p_corr_ep", std::string(p1[k]) + "_theta", std::string(p1[k]) + "_phi",
                std::string(p2[k]) + "_p_corr_em", std::string(p2[k]) + "_theta", std::string(p2[k]) + "_phi",
                "pip_p_corr_pip", "pip_theta", "pip_phi", "pim_p_corr_pim", "pim_theta", "pim_phi"};
            std::vector<Float_t> x(nm.size());
            for (size_t i = 0; i < nm.size(); ++i) { t->SetBranchStatus(nm[i].c_str(), 1); t->SetBranchAddress(nm[i].c_str(), &x[i]); }
            const Long64_t n = t->GetEntries();
            for (Long64_t e = 0; e < n; ++e) {
                t->GetEntry(e);
                if (!(x[0] > -0.5)) continue;
                const float vz = x[2];
                const bool z500 = vz > -500.f, pt3 = (x[1] == 8192.f);
                TLorentzVector sum; bool cor = true;
                const double ms[4] = {me, me, mpi, mpi};
                for (int q = 0; q < 4; ++q) {
                    const double p = x[4 + 3 * q], th = x[5 + 3 * q] * TMath::DegToRad(), ph = x[6 + 3 * q] * TMath::DegToRad();
                    if (!std::isfinite(p)) { cor = false; break; }
                    TLorentzVector v; v.SetXYZM(p * std::sin(th) * std::cos(ph), p * std::sin(th) * std::sin(ph), p * std::cos(th), ms[q]);
                    sum += v;
                }
                c[k][0]++;
                if (z500) c[k][1]++;
                if (z500 && pt3) c[k][2]++;
                if (z500 && pt3 && cor) c[k][3]++;
                if (z500 && pt3 && cor && vz > -160.f) c[k][4]++;
                if (!z500) { c[k][5]++; if (vz == -1000.f) c[k][6]++; if (x[3] == 0.f) c[k][7]++; }
                if (cor) {
                    const double m4 = 1e-3 * sum.M();
                    if (m4 >= pkLo && m4 < pkHi) { c[k][8]++; if (z500) c[k][9]++; }
                    if (m4 >= bpLo && m4 < bpHi) { c[k][10]++; if (z500) c[k][11]++; }
                }
            }
        }
        f->Close(); delete f;
        fprintf(stderr, "file %s done\n", files[fi].c_str());
    }
    for (int k = 0; k < 3; ++k) {
        printf("VZ500 %d", k);
        for (int i = 0; i < 12; ++i) printf(" %lld", c[k][i]);
        printf("\n");
    }
}
