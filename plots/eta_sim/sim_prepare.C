// ============================================================================
// sim_prepare.C  --  compact simulation file for the eta study in M(pi+pi-e+e-)
// ============================================================================
// One pass over scan_nt of the FAT outputs of the SMASH simulation.
//   part 0 (default, the first study): the six std and lep outputs
//     /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/output_{std,lep}_{pippimepem,pippimepep,pippimemem}.root
//   writes
//     /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/sim_compact.root
//          trees std_us, std_pp, std_mm, lep_us, lep_pp, lep_mm (columns = simeng::Row)
//     /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/sim_prepare.log
//   part 1: the three eta-set outputs (sim_set = 3)
//     /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/output_eta_{pippimepem,pippimepep,pippimemem}.root
//   writes
//     /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/sim_compact_eta.root
//          trees eta_us, eta_pp, eta_mm (the same columns)
//     /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/sim_prepare_eta.log
//   outPath (optional): write the compact file elsewhere (log = outPath + ".log"); used to
//     check that part 0 still reproduces sim_compact.root column by column
//
// Every row of scan_nt is kept (the FAT selection: isBest >= 0, z > -500,
// trigger bit 8192); vz > -160 and the other cuts are applied by the engine.
// Per row it adds to the stored values, with the code of ../scan_prepare.C so
// that names and definitions equal the data compact file:
//   - event group (maximal run of adjacent rows with equal (event, eVertReco_z);
//     runnumber is 0 and file_id -1 in simulation): grp (unique per tree), fold
//     (grp parity), nrow and the number of distinct tracks per slot
//   - crossing points of the two leptons and of the two pions
//   - the photon whose M(e+e-gamma) is closest to 0.140 GeV/c2 (m_eeg, m5, mm5)
//   - T_pi0 and the kinematic limit M4max for the measured M(pi+pi-)
//   - pair quantities of the pi+pi- and e+e- subsystems (simeng::ComputePairKin)
//   - the truth: classes, template, weights (double copies), ownership, flags,
//     ghost bits, GEANT ids and parents, true masses, true opening angle,
//     true momenta, polar angles and creation vertices
//
// Run from /home/przygoda/HADES/FAT/FAT (ACLiC, build directory outside the repository):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_engine", kTRUE)' 'plots/eta_sim/sim_prepare.C+'
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_eta", kTRUE)' 'plots/eta_sim/sim_prepare.C+(1)'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TVector3.h"

#include "SimEngine.h"

namespace sprep {

#define SPREP_TRACK_COLS(X, T)                                                   \
    X(T##_p_rec) X(T##_p_cor) X(T##_theta) X(T##_phi) X(T##_system) X(T##_r)     \
    X(T##_z) X(T##_rkchi2) X(T##_isUsedVertex) X(T##_oa_lept) X(T##_oa_hadr)     \
    X(T##_sim_id) X(T##_sim_parentid) X(T##_sim_grandparentid)                   \
    X(T##_sim_processid) X(T##_sim_p) X(T##_sim_px) X(T##_sim_py) X(T##_sim_pz)  \
    X(T##_sim_vertexx) X(T##_sim_vertexy) X(T##_sim_vertexz) X(T##_sim_isghost)

#define SPREP_IN_COLS(X)                                                         \
    X(pat_entry) X(event) X(runnumber) X(trigbit) X(isBest) X(hyp_chi2)          \
    X(eVertReco_x) X(eVertReco_y) X(eVertReco_z) X(eVertReco_chi2)               \
    X(neutr_mult) X(fwdet_mult) X(evtStart) X(cor_ok)                            \
    SPREP_TRACK_COLS(X, ep) SPREP_TRACK_COLS(X, em)                              \
    SPREP_TRACK_COLS(X, pip) SPREP_TRACK_COLS(X, pim)                            \
    X(ep_rich_ind) X(em_rich_ind) X(ep_richmatchquality) X(em_richmatchquality)  \
    X(ep_ringchi2) X(em_ringchi2) X(ep_kIsLepton) X(em_kIsLepton)                \
    X(ep_rich_sim_geanttrack) X(em_rich_sim_geanttrack)                          \
    X(pip_isring) X(pim_isring) X(pip_mass2) X(pim_mass2)                        \
    X(pip_dedx_mdc) X(pim_dedx_mdc)                                              \
    X(oa) X(m_pippimepem_cor) X(m_ee_cor) X(m_pippim_cor)                        \
    X(oa_pippim_epem_lab_cor) X(oa_pippim_epem_eta_rest_cor)                     \
    X(mm_pippimepem_cor) X(mm_pippim_cor) X(mm_epem_cor) X(pt4_cor) X(y4_cms_cor) \
    X(n_gam) X(g1_e) X(g1_theta) X(g1_phi) X(g2_e) X(g2_theta) X(g2_phi)         \
    X(g3_e) X(g3_theta) X(g3_phi) X(fw1_p) X(fw1_theta)                          \
    X(sim_set) X(sim_cls) X(sim_lcls) X(sim_pcls) X(sim_tmpl) X(sim_w)           \
    X(sim_w_db) X(sim_own_A) X(sim_own_B) X(sim_flags)                           \
    X(m4_true) X(mee_true) X(mpipi_true)

struct In {
#define X(n) Float_t n;
    SPREP_IN_COLS(X)
#undef X
};

struct Trk { float p, th, ph; };
inline bool Same(const Trk& a, const Trk& b) { return a.p == b.p && a.th == b.th && a.ph == b.ph; }
inline void AddDistinct(std::vector<Trk>& v, const Trk& t) {
    for (const auto& x : v) if (Same(x, t)) return;
    v.push_back(t);
}
inline int I(float x) { return (int)std::lround(x); }

struct Stat {
    Long64_t rows = 0, groups = 0, best1 = 0, grpBestNot1 = 0, corbad = 0, nanM4 = 0, nanW = 0;
    Long64_t rowsVz = 0, ownAvz = 0, withGam = 0, parallelEE = 0, parallelPP = 0;
    double sumWvz = 0, sumWdbVz = 0;          // owned in mode A, vz > -160 (double)
    double maxd[3] = {0, 0, 0};               // |recomputed - stored| of m4, mee, mpipi
    Long64_t nrowHist[6] = {0, 0, 0, 0, 0, 0};
    Long64_t multiPip = 0, multiPim = 0, multiL1 = 0, multiL2 = 0;
};

}  // namespace sprep

void sim_prepare(int part = 0, const char* outPath = nullptr) {
    using namespace sprep;
    using namespace simeng;
    TStopwatch swAll;
    if (part != 0 && part != 1) { printf("sim_prepare: part must be 0 (std + lep) or 1 (eta)\n"); return; }
    const char* inName[kNTreeAll] = {"output_std_pippimepem.root", "output_std_pippimepep.root",
                                     "output_std_pippimemem.root", "output_lep_pippimepem.root",
                                     "output_lep_pippimepep.root", "output_lep_pippimemem.root",
                                     "output_eta_pippimepem.root", "output_eta_pippimepep.root",
                                     "output_eta_pippimemem.root"};
    // expected numbers (validated productions of 2026-11-12: tables/engine_validate_summary.txt for std and lep,
    // tables/eta_validate.txt (sim_eta_checks.C) for eta); every group has exactly one isBest == 1 row
    const Long64_t expRows[kNTreeAll] = {230469, 1133, 4791, 7441490, 205527, 495711, 3366025, 13881, 53628};
    const Long64_t expGrp[kNTreeAll]  = {213044, 1053, 4397, 6856223, 196123, 466029, 2731582, 12233, 45375};
    const Long64_t expVz[kNTreeAll]   = {230115, 1132, 4777, 7432170, 205288, 495101, 3357676, 13856, 53510};
    const int kFirst = part == 0 ? 0 : kNTree, kLast = part == 0 ? kNTree : kNTreeAll;
    const std::string outName = outPath ? std::string(outPath) : std::string(part == 0 ? kCompact : kCompactEta);
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const TLorentzVector ini = Initial();

    TFile* fout = TFile::Open(outName.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) { printf("cannot create %s\n", outName.c_str()); return; }
    fout->SetCompressionSettings(404);   // LZ4: fast to read

    Stat st[kNTreeAll];
    bool ok = true;

    for (int k = kFirst; k < kLast; ++k) {
        TStopwatch sw;
        const int setId = k / 3 + 1;   // 1 std, 2 lep, 3 eta
        const std::string fn = std::string(kProdDir) + "/" + inName[k];
        TFile* fin = TFile::Open(fn.c_str(), "READ");
        if (!fin || fin->IsZombie()) { printf("cannot open %s\n", fn.c_str()); return; }
        TTree* tin = (TTree*)fin->Get("scan_nt");
        if (!tin) { printf("scan_nt missing in %s\n", fn.c_str()); return; }
        In in;
        tin->SetBranchStatus("*", 0);
#define X(n) tin->SetBranchStatus(#n, 1); tin->SetBranchAddress(#n, &in.n);
        SPREP_IN_COLS(X)
#undef X
        tin->SetCacheSize(300 * 1024 * 1024);

        fout->cd();
        TTree* tout = new TTree(TreeName(k), Form("compact simulation rows, %s", inName[k]));
        Row row;
        BindRow(tout, row, true);
        tout->SetAutoSave(0);

        std::vector<Row> buf;
        std::vector<Trk> tPip, tPim, tL1, tL2, tLall;
        int nBestInGroup = 0;
        Long64_t grp = -1;
        float keyE = -1, keyZ = -1e9;
        Stat& S = st[k];

        auto flush = [&]() {
            if (buf.empty()) return;
            const int nrow = (int)buf.size();
            for (auto& r : buf) {
                r.nrow = nrow;
                r.n_pip = (int)tPip.size(); r.n_pim = (int)tPim.size();
                r.n_lep1 = (int)tL1.size(); r.n_lep2 = (int)tL2.size(); r.n_lep = (int)tLall.size();
                row = r;
                tout->Fill();
            }
            S.nrowHist[nrow >= 6 ? 5 : nrow - 1]++;
            if (tPip.size() > 1) S.multiPip++;
            if (tPim.size() > 1) S.multiPim++;
            if (tL1.size() > 1) S.multiL1++;
            if (tL2.size() > 1) S.multiL2++;
            if (nBestInGroup != 1) S.grpBestNot1++;
            buf.clear(); tPip.clear(); tPim.clear(); tL1.clear(); tL2.clear(); tLall.clear();
            nBestInGroup = 0;
        };

        const Long64_t n = tin->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            tin->GetEntry(i);
            if (i == 0 || in.event != keyE || in.eVertReco_z != keyZ) {
                flush();
                ++grp;
                keyE = in.event; keyZ = in.eVertReco_z;
            }
            Row r;
            // ---- keys ----
            r.set = setId;
            r.event = I(in.event);
            r.pat_entry = I(in.pat_entry);
            r.trigbit = I(in.trigbit);
            r.grp = (int)grp;
            r.fold = (int)(grp & 1);
            r.isBest = I(in.isBest);
            r.cor_ok = I(in.cor_ok);
            r.nrow = r.n_pip = r.n_pim = r.n_lep1 = r.n_lep2 = r.n_lep = 0;
            if (r.isBest == 1) { nBestInGroup++; S.best1++; }
            if (!r.cor_ok) S.corbad++;

            // ---- stored corrected values (identical to the data compact file) ----
            r.m4 = in.m_pippimepem_cor; r.mee = in.m_ee_cor; r.mpipi = in.m_pippim_cor;
            r.oa = in.oa; r.oa_lab = in.oa_pippim_epem_lab_cor; r.oa_eta = in.oa_pippim_epem_eta_rest_cor;
            r.mm4 = in.mm_pippimepem_cor; r.mm_pipi = in.mm_pippim_cor; r.mm_ee = in.mm_epem_cor;
            r.pt4 = in.pt4_cor; r.y4 = in.y4_cms_cor;
            if (r.m4 != r.m4) S.nanM4++;

            // ---- four-vectors: corrected momentum where it is a number, reconstructed otherwise (data rule) ----
            const double pEp  = std::isnan(in.ep_p_cor)  ? in.ep_p_rec  : in.ep_p_cor;
            const double pEm  = std::isnan(in.em_p_cor)  ? in.em_p_rec  : in.em_p_cor;
            const double pPip = std::isnan(in.pip_p_cor) ? in.pip_p_rec : in.pip_p_cor;
            const double pPim = std::isnan(in.pim_p_cor) ? in.pim_p_rec : in.pim_p_cor;
            const TLorentzVector vEp  = P4(pEp,  in.ep_theta,  in.ep_phi,  kMeMeV);
            const TLorentzVector vEm  = P4(pEm,  in.em_theta,  in.em_phi,  kMeMeV);
            const TLorentzVector vPip = P4(pPip, in.pip_theta, in.pip_phi, kMpiMeV);
            const TLorentzVector vPim = P4(pPim, in.pim_theta, in.pim_phi, kMpiMeV);
            const TLorentzVector vEE = vEp + vEm, vPP = vPip + vPim, v4 = vEE + vPP;
            if (r.cor_ok) {
                const double d0 = std::fabs(v4.M() / 1000.0 - r.m4), d1 = std::fabs(vEE.M() / 1000.0 - r.mee),
                             d2 = std::fabs(vPP.M() / 1000.0 - r.mpipi);
                if (d0 > S.maxd[0]) S.maxd[0] = d0;
                if (d1 > S.maxd[1]) S.maxd[1] = d1;
                if (d2 > S.maxd[2]) S.maxd[2] = d2;
            }
            // pair quantities (corrected momenta; NaN when cor_ok == 0, as every corrected variable)
            {
                const PairKin pk = ComputePairKin(in.ep_p_cor, in.ep_theta, in.ep_phi, in.em_p_cor, in.em_theta,
                                                  in.em_phi, in.pip_p_cor, in.pip_theta, in.pip_phi, in.pim_p_cor,
                                                  in.pim_theta, in.pim_phi);
                r.p_pipi = pk.p_pipi; r.th_pipi = pk.th_pipi; r.pt_pipi = pk.pt_pipi; r.y_pipi = pk.y_pipi;
                r.p_ee = pk.p_ee; r.th_ee = pk.th_ee; r.pt_ee = pk.pt_ee; r.y_ee = pk.y_ee;
            }

            // ---- vertex, tracks ----
            r.vx = in.eVertReco_x; r.vy = in.eVertReco_y; r.vz = in.eVertReco_z; r.vchi2 = in.eVertReco_chi2;
            r.ep_p = in.ep_p_cor; r.em_p = in.em_p_cor; r.pip_p = in.pip_p_cor; r.pim_p = in.pim_p_cor;
            r.ep_p_rec = in.ep_p_rec; r.em_p_rec = in.em_p_rec; r.pip_p_rec = in.pip_p_rec; r.pim_p_rec = in.pim_p_rec;
            r.ep_theta = in.ep_theta; r.em_theta = in.em_theta; r.pip_theta = in.pip_theta; r.pim_theta = in.pim_theta;
            r.ep_phi = in.ep_phi; r.em_phi = in.em_phi; r.pip_phi = in.pip_phi; r.pim_phi = in.pim_phi;
            r.ep_z = in.ep_z; r.em_z = in.em_z; r.pip_z = in.pip_z; r.pim_z = in.pim_z;
            r.ep_r = in.ep_r; r.em_r = in.em_r; r.pip_r = in.pip_r; r.pim_r = in.pim_r;
            r.ep_rkchi2 = in.ep_rkchi2; r.em_rkchi2 = in.em_rkchi2; r.pip_rkchi2 = in.pip_rkchi2; r.pim_rkchi2 = in.pim_rkchi2;
            r.ep_system = I(in.ep_system); r.em_system = I(in.em_system);
            r.pip_system = I(in.pip_system); r.pim_system = I(in.pim_system);
            r.ep_isUsedVertex = I(in.ep_isUsedVertex); r.em_isUsedVertex = I(in.em_isUsedVertex);
            r.pip_isUsedVertex = I(in.pip_isUsedVertex); r.pim_isUsedVertex = I(in.pim_isUsedVertex);
            r.ep_kIsLepton = I(in.ep_kIsLepton); r.em_kIsLepton = I(in.em_kIsLepton);
            r.ep_richmatchquality = in.ep_richmatchquality; r.em_richmatchquality = in.em_richmatchquality;
            r.ep_ringchi2 = in.ep_ringchi2; r.em_ringchi2 = in.em_ringchi2;
            r.ep_oa_lept = in.ep_oa_lept; r.em_oa_lept = in.em_oa_lept;
            r.ep_oa_hadr = in.ep_oa_hadr; r.em_oa_hadr = in.em_oa_hadr;
            r.shared_ring = (in.ep_rich_ind >= 0 && in.ep_rich_ind == in.em_rich_ind) ? 1 : 0;
            r.pip_isring = I(in.pip_isring); r.pim_isring = I(in.pim_isring);
            r.pip_mass2 = in.pip_mass2; r.pim_mass2 = in.pim_mass2;
            r.pip_dedx_mdc = in.pip_dedx_mdc; r.pim_dedx_mdc = in.pim_dedx_mdc;
            r.evtStart = I(in.evtStart);
            r.hyp_chi2 = in.hyp_chi2;
            r.neutr_mult = I(in.neutr_mult);
            r.fwdet_mult = I(in.fwdet_mult);
            r.fw1_p = in.fw1_p; r.fw1_theta = in.fw1_theta;

            // ---- crossing points ----
            const CrossPt cE = Crossing(in.ep_r, in.ep_z, in.ep_theta, in.ep_phi, in.em_r, in.em_z, in.em_theta,
                                        in.em_phi, in.eVertReco_x, in.eVertReco_y, in.eVertReco_z);
            const CrossPt cP = Crossing(in.pip_r, in.pip_z, in.pip_theta, in.pip_phi, in.pim_r, in.pim_z,
                                        in.pim_theta, in.pim_phi, in.eVertReco_x, in.eVertReco_y, in.eVertReco_z);
            r.dz_ee = cE.dz; r.r_ee = cE.r; r.rv_ee = cE.rv; r.doca_ee = cE.doca;
            r.dz_pipi = cP.dz; r.r_pipi = cP.r; r.rv_pipi = cP.rv; r.doca_pipi = cP.doca;
            if (cE.r < 0) S.parallelEE++;
            if (cP.r < 0) S.parallelPP++;

            // ---- photon with M(e+e-gamma) closest to 0.140 (as ../scan_prepare.C) ----
            r.n_gam = I(in.n_gam);
            r.g_idx = -1; r.g_e = r.g_theta = r.g_phi = -1.f;
            r.m_eeg = r.m5 = r.mm5 = -1.f;
            {
                const float ge[3] = {in.g1_e, in.g2_e, in.g3_e};
                const float gt[3] = {in.g1_theta, in.g2_theta, in.g3_theta};
                const float gp[3] = {in.g1_phi, in.g2_phi, in.g3_phi};
                double best = 1e9;
                TLorentzVector vBest;
                for (int j = 0; j < 3; ++j) {
                    if (!(ge[j] > 0)) continue;
                    const TLorentzVector vg = P4(ge[j], gt[j], gp[j], 0.0);
                    const double m = (vEE + vg).M() / 1000.0;
                    if (std::fabs(m - 0.140) < best) {
                        best = std::fabs(m - 0.140);
                        r.g_idx = j + 1; r.g_e = ge[j]; r.g_theta = gt[j]; r.g_phi = gp[j];
                        r.m_eeg = (float)m;
                        vBest = vg;
                    }
                }
                if (r.g_idx > 0) {
                    S.withGam++;
                    r.m5 = r.cor_ok ? (float)((v4 + vBest).M() / 1000.0) : nan;
                    r.mm5 = r.cor_ok ? (float)((ini - v4 - vBest).M() / 1000.0) : nan;
                }
            }

            // ---- derived (as data) ----
            r.T_pi0 = (float)scan::TPi0(r.mpipi);
            r.m4max = std::isnan(r.mpipi) ? nan : (float)scan::M4max(r.mpipi);

            // ---- truth ----
            r.cls = I(in.sim_cls); r.lcls = I(in.sim_lcls); r.pcls = I(in.sim_pcls); r.tmpl = I(in.sim_tmpl);
            r.ownA = I(in.sim_own_A); r.ownB = I(in.sim_own_B); r.flags = I(in.sim_flags);
            r.w = (double)in.sim_w; r.w_db = (double)in.sim_w_db;
            if (r.w != r.w || r.w_db != r.w_db) S.nanW++;
            r.ghost = (I(in.ep_sim_isghost) == 1 ? 1 : 0) | (I(in.em_sim_isghost) == 1 ? 2 : 0) |
                      (I(in.pip_sim_isghost) == 1 ? 4 : 0) | (I(in.pim_sim_isghost) == 1 ? 8 : 0);
            r.ep_sim_id = I(in.ep_sim_id); r.em_sim_id = I(in.em_sim_id);
            r.pip_sim_id = I(in.pip_sim_id); r.pim_sim_id = I(in.pim_sim_id);
            r.ep_sim_parentid = I(in.ep_sim_parentid); r.em_sim_parentid = I(in.em_sim_parentid);
            r.pip_sim_parentid = I(in.pip_sim_parentid); r.pim_sim_parentid = I(in.pim_sim_parentid);
            r.ep_sim_grandparentid = I(in.ep_sim_grandparentid); r.em_sim_grandparentid = I(in.em_sim_grandparentid);
            r.pip_sim_grandparentid = I(in.pip_sim_grandparentid); r.pim_sim_grandparentid = I(in.pim_sim_grandparentid);
            r.ep_sim_processid = I(in.ep_sim_processid); r.em_sim_processid = I(in.em_sim_processid);
            r.pip_sim_processid = I(in.pip_sim_processid); r.pim_sim_processid = I(in.pim_sim_processid);
            r.ep_rich_sim_geanttrack = I(in.ep_rich_sim_geanttrack);
            r.em_rich_sim_geanttrack = I(in.em_rich_sim_geanttrack);
            r.m4_true = in.m4_true; r.mee_true = in.mee_true; r.mpipi_true = in.mpipi_true;
            {
                const float px[4] = {in.ep_sim_px, in.em_sim_px, in.pip_sim_px, in.pim_sim_px};
                const float py[4] = {in.ep_sim_py, in.em_sim_py, in.pip_sim_py, in.pim_sim_py};
                const float pz[4] = {in.ep_sim_pz, in.em_sim_pz, in.pip_sim_pz, in.pim_sim_pz};
                const float vxx[4] = {in.ep_sim_vertexx, in.em_sim_vertexx, in.pip_sim_vertexx, in.pim_sim_vertexx};
                const float vyy[4] = {in.ep_sim_vertexy, in.em_sim_vertexy, in.pip_sim_vertexy, in.pim_sim_vertexy};
                const float vzz[4] = {in.ep_sim_vertexz, in.em_sim_vertexz, in.pip_sim_vertexz, in.pim_sim_vertexz};
                const float sp[4] = {in.ep_sim_p, in.em_sim_p, in.pip_sim_p, in.pim_sim_p};
                float th[4], vr[4];
                TVector3 pv[4];
                for (int j = 0; j < 4; ++j) {
                    pv[j].SetXYZ(px[j], py[j], pz[j]);
                    th[j] = pv[j].Mag() > 0 ? (float)(pv[j].Theta() * kR2D) : -1.f;
                    vr[j] = (float)std::hypot(vxx[j], vyy[j]);
                }
                r.ep_sim_p = sp[0]; r.em_sim_p = sp[1]; r.pip_sim_p = sp[2]; r.pim_sim_p = sp[3];
                r.ep_sim_theta = th[0]; r.em_sim_theta = th[1]; r.pip_sim_theta = th[2]; r.pim_sim_theta = th[3];
                r.ep_sim_vz = vzz[0]; r.em_sim_vz = vzz[1]; r.pip_sim_vz = vzz[2]; r.pim_sim_vz = vzz[3];
                r.ep_sim_vr = vr[0]; r.em_sim_vr = vr[1]; r.pip_sim_vr = vr[2]; r.pim_sim_vr = vr[3];
                r.oa_true = (pv[0].Mag() > 0 && pv[1].Mag() > 0) ? (float)(pv[0].Angle(pv[1]) * kR2D) : -1.f;
            }
            if (r.vz > -160.f) {
                S.rowsVz++;
                if (r.ownA == 1) { S.ownAvz++; S.sumWvz += r.w; S.sumWdbVz += r.w_db; }
            }

            // ---- distinct tracks of the group (as data) ----
            AddDistinct(tPip, Trk{in.pip_p_rec, in.pip_theta, in.pip_phi});
            AddDistinct(tPim, Trk{in.pim_p_rec, in.pim_theta, in.pim_phi});
            const Trk c{in.ep_p_rec, in.ep_theta, in.ep_phi}, d{in.em_p_rec, in.em_theta, in.em_phi};
            AddDistinct(tL1, c); AddDistinct(tL2, d); AddDistinct(tLall, c); AddDistinct(tLall, d);

            buf.push_back(r);
            S.rows++;
        }
        flush();
        S.groups = grp + 1;
        fout->cd();
        tout->Write("", TObject::kOverwrite);
        fin->Close();
        delete fin;
        printf("[%s] %lld rows, %lld groups, %.1f s\n", TreeName(k), S.rows, S.groups, sw.RealTime());
        if (S.rows != expRows[k] || S.groups != expGrp[k] || S.best1 != expGrp[k] || S.grpBestNot1 != 0 ||
            S.rowsVz != expVz[k] || S.nanW != 0)
            ok = false;
    }
    fout->Close();
    delete fout;

    // ---- read-back: compact file against the source, row by row ----
    Long64_t rbRows[kNTreeAll] = {0}, rbGroups[kNTreeAll] = {0}, rbBest[kNTreeAll] = {0}, rbBad[kNTreeAll] = {0};
    double rbSumW[kNTreeAll] = {0};
    {
        TFile* fc = TFile::Open(outName.c_str(), "READ");
        for (int k = kFirst; k < kLast && fc && !fc->IsZombie(); ++k) {
            TTree* tc = (TTree*)fc->Get(TreeName(k));
            TFile* fin = TFile::Open((std::string(kProdDir) + "/" + inName[k]).c_str(), "READ");
            TTree* tin = fin ? (TTree*)fin->Get("scan_nt") : nullptr;
            if (!tc || !tin) { ok = false; continue; }
            Int_t cGrp = 0, cBest = 0, cPat = 0, cOwn = 0, cCls = 0;
            Float_t cM4 = 0, cVz = 0, sM4 = 0, sPat = 0, sW = 0, sCls = 0;
            Double_t cW = 0;
            tc->SetBranchStatus("*", 0);
            for (const char* b : {"grp", "isBest", "pat_entry", "ownA", "cls", "m4", "vz", "w"}) tc->SetBranchStatus(b, 1);
            tc->SetBranchAddress("grp", &cGrp); tc->SetBranchAddress("isBest", &cBest);
            tc->SetBranchAddress("pat_entry", &cPat); tc->SetBranchAddress("ownA", &cOwn);
            tc->SetBranchAddress("cls", &cCls); tc->SetBranchAddress("m4", &cM4);
            tc->SetBranchAddress("vz", &cVz); tc->SetBranchAddress("w", &cW);
            tin->SetBranchStatus("*", 0);
            for (const char* b : {"m_pippimepem_cor", "pat_entry", "sim_w", "sim_cls"}) tin->SetBranchStatus(b, 1);
            tin->SetBranchAddress("m_pippimepem_cor", &sM4); tin->SetBranchAddress("pat_entry", &sPat);
            tin->SetBranchAddress("sim_w", &sW); tin->SetBranchAddress("sim_cls", &sCls);
            const Long64_t n = tc->GetEntries();
            rbRows[k] = n;
            if (n != tin->GetEntries()) ok = false;
            Int_t lastGrp = -1;
            for (Long64_t i = 0; i < n; ++i) {
                tc->GetEntry(i); tin->GetEntry(i);
                if (cGrp != lastGrp) { rbGroups[k]++; lastGrp = cGrp; }
                if (cBest == 1) rbBest[k]++;
                if (cPat != (int)std::lround(sPat) || !(cM4 == sM4) || cW != (double)sW || cCls != (int)std::lround(sCls)) rbBad[k]++;
                if (cVz > -160.f && cOwn == 1) rbSumW[k] += cW;
            }
            if (rbRows[k] != expRows[k] || rbGroups[k] != expGrp[k] || rbBest[k] != expGrp[k] || rbBad[k] != 0) ok = false;
            fin->Close();
            delete fin;
        }
        if (fc) { fc->Close(); delete fc; }
    }

    // ---- verification report ----
    const std::string logName = outPath ? outName + ".log"
                                        : std::string(kScanDir) + (part == 0 ? "/sim_prepare.log" : "/sim_prepare_eta.log");
    FILE* lf = fopen(logName.c_str(), "w");
    for (int pass = 0; pass < 2; ++pass) {
        FILE* o = pass == 0 ? stdout : lf;
        if (!o) continue;
        fprintf(o, "sim_prepare: verification of %s\n", outName.c_str());
        for (int k = kFirst; k < kLast; ++k) {
            const Stat& S = st[k];
            fprintf(o, "\n[%s] %s\n", TreeName(k), inName[k]);
            fprintf(o, "  rows %lld (expected %lld)  groups %lld (expected %lld)  isBest==1 rows %lld\n", S.rows,
                    expRows[k], S.groups, expGrp[k], S.best1);
            fprintf(o, "  groups without exactly one isBest==1 row: %lld ; cor_ok==0 rows %lld ; NaN m4 %lld ; NaN w %lld\n",
                    S.grpBestNot1, S.corbad, S.nanM4, S.nanW);
            fprintf(o, "  vz > -160: rows %lld (expected %lld), owned in mode A %lld, sum w %.8g, sum w_db %.8g\n",
                    S.rowsVz, expVz[k], S.ownAvz, S.sumWvz, S.sumWdbVz);
            fprintf(o, "  read-back: rows %lld, groups %lld, isBest==1 %lld, rows differing from the source %lld, sum w (owned A, vz > -160) %.8g\n",
                    rbRows[k], rbGroups[k], rbBest[k], rbBad[k], rbSumW[k]);
            fprintf(o, "  max |recomputed - stored| m4 %.2g, mee %.2g, mpipi %.2g GeV/c2\n", S.maxd[0], S.maxd[1], S.maxd[2]);
            fprintf(o, "  groups with 1,2,3,4,5,>=6 rows: %lld %lld %lld %lld %lld %lld\n", S.nrowHist[0], S.nrowHist[1],
                    S.nrowHist[2], S.nrowHist[3], S.nrowHist[4], S.nrowHist[5]);
            fprintf(o, "  groups with more than one distinct pi+ %lld, pi- %lld, lepton1 %lld, lepton2 %lld\n", S.multiPip,
                    S.multiPim, S.multiL1, S.multiL2);
            fprintf(o, "  rows with a photon %lld ; parallel lines: leptons %lld, pions %lld\n", S.withGam, S.parallelEE,
                    S.parallelPP);
        }
        fprintf(o, "\nresult: %s\n", ok ? "row, group and weight checks as expected" : "MISMATCH (see above)");
    }
    if (lf) fclose(lf);
    printf("sim_prepare: total %.1f s, log %s\n", swAll.RealTime(), logName.c_str());
}
