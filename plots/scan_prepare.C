// ============================================================================
// scan_prepare.C  --  compact data file for the eta cut scan in M(pi+pi-e+e-)
// ============================================================================
// One pass over scan_nt of the three hypotheses (e+e-, e+e+, e-e-) of the
// production pippimepem_isBestAll_T4530. Writes
//   <prod>/scan/scan_compact.root   trees us, pp, mm (columns = scan::Row)
//   <prod>/scan/scan_prepare.log    verification numbers
//
// Per row it adds to the stored values:
//   - event group (maximal run of adjacent rows with equal file_id, runnumber,
//     event): grp, fold, nrow and the number of distinct tracks per slot
//   - the fallback set *_fb (reconstructed pi- momentum when cor_ok == 0)
//   - crossing points of the two leptons and of the two pions
//   - the photon whose M(e+e-gamma) is closest to 0.140 GeV/c2
//   - T_pi0 and the kinematic limit M4max for the measured M(pi+pi-)
//
// Run (ACLiC, build directory outside the repository):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan/aclic_engine", kTRUE)' 'plots/scan_prepare.C+'
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

#include "ScanEngine.h"

namespace prep {

// ---- input columns of scan_nt (all Float_t) ---------------------------------
#define PREP_TRACK_COLS(X, T)                                                    \
    X(T##_p_rec) X(T##_p_cor) X(T##_theta) X(T##_phi) X(T##_system) X(T##_r)     \
    X(T##_z) X(T##_rkchi2) X(T##_isUsedVertex) X(T##_oa_lept) X(T##_oa_hadr)

#define PREP_IN_COLS(X)                                                          \
    X(pat_entry) X(file_id) X(runnumber) X(event) X(isBest) X(hyp_chi2)          \
    X(eVertReco_x) X(eVertReco_y) X(eVertReco_z) X(eVertReco_chi2)               \
    X(neutr_mult) X(fwdet_mult) X(evtStart) X(cor_ok)                            \
    PREP_TRACK_COLS(X, ep) PREP_TRACK_COLS(X, em)                                \
    PREP_TRACK_COLS(X, pip) PREP_TRACK_COLS(X, pim)                              \
    X(ep_rich_ind) X(em_rich_ind) X(ep_richmatchquality) X(em_richmatchquality)  \
    X(ep_ringchi2) X(em_ringchi2) X(ep_kIsLepton) X(em_kIsLepton)                \
    X(pip_isring) X(pim_isring) X(pip_mass2) X(pim_mass2)                        \
    X(pip_dedx_mdc) X(pim_dedx_mdc)                                              \
    X(oa) X(m_pippimepem_cor) X(m_ee_cor) X(m_pippim_cor)                        \
    X(oa_pippim_epem_lab_cor) X(oa_pippim_epem_eta_rest_cor)                     \
    X(mm_pippimepem_cor) X(mm_pippim_cor) X(mm_epem_cor) X(pt4_cor) X(y4_cms_cor) \
    X(n_gam) X(g1_e) X(g1_theta) X(g1_phi) X(g2_e) X(g2_theta) X(g2_phi)         \
    X(g3_e) X(g3_theta) X(g3_phi) X(fw1_p) X(fw1_theta)

struct In {
#define X(n) Float_t n;
    PREP_IN_COLS(X)
#undef X
};

const double kD2R = 1.74532925199432955e-02;
const double kR2D = 57.2957795130823229;

TLorentzVector P4(double p, double theta, double phi, double mMeV) {
    const double th = theta * kD2R, ph = phi * kD2R;
    TLorentzVector v;
    v.SetVectM(TVector3(p * std::sin(th) * std::cos(ph), p * std::sin(th) * std::sin(ph),
                        p * std::cos(th)), mMeV);
    return v;
}

// Point of closest approach of two straight lines. Track line as in HYDRA
// (HParticleTool::calcSegVector): base = (r cos(phi+90), r sin(phi+90), z),
// dir = (sin(theta) cos(phi), sin(theta) sin(phi), cos(theta)); r is signed.
struct CrossPt { float dz, r, rv, doca; };
CrossPt Crossing(double r1, double z1, double th1, double ph1,
               double r2, double z2, double th2, double ph2,
               double vx, double vy, double vz) {
    const double t1 = th1 * kD2R, p1 = ph1 * kD2R, t2 = th2 * kD2R, p2 = ph2 * kD2R;
    const TVector3 b1(-r1 * std::sin(p1), r1 * std::cos(p1), z1);
    const TVector3 b2(-r2 * std::sin(p2), r2 * std::cos(p2), z2);
    const TVector3 d1(std::sin(t1) * std::cos(p1), std::sin(t1) * std::sin(p1), std::cos(t1));
    const TVector3 d2(std::sin(t2) * std::cos(p2), std::sin(t2) * std::sin(p2), std::cos(t2));
    const TVector3 w = b1 - b2;
    const double b = d1.Dot(d2), d = d1.Dot(w), e = d2.Dot(w);
    const double den = 1.0 - b * b;
    CrossPt c;
    if (den < 1e-12) {                       // parallel lines: no crossing point
        c.dz = -9999.f; c.r = -1.f; c.rv = -1.f;
        c.doca = (float)(w - d * d1).Mag();
        return c;
    }
    const double s = (b * e - d) / den, t = (e - b * d) / den;
    const TVector3 q1 = b1 + s * d1, q2 = b2 + t * d2;
    const TVector3 q = 0.5 * (q1 + q2);
    c.dz = (float)(q.Z() - vz);
    c.r = (float)std::hypot(q.X(), q.Y());
    c.rv = (float)std::hypot(q.X() - vx, q.Y() - vy);
    c.doca = (float)(q1 - q2).Mag();
    return c;
}

// distance of a track line from a point (used only for the convention check)
double DistToPoint(double r, double z, double th, double ph, double x, double y, double zz, int conv) {
    const double t = th * kD2R, p = ph * kD2R;
    const double sg = conv == 0 ? 1.0 : -1.0;     // conv 1 = opposite sign of r
    const TVector3 b(-sg * r * std::sin(p), sg * r * std::cos(p), z);
    const TVector3 d(std::sin(t) * std::cos(p), std::sin(t) * std::sin(p), std::cos(t));
    const TVector3 w = TVector3(x, y, zz) - b;
    return (w - w.Dot(d) * d).Mag();
}

struct Trk { float p, th, ph; };
bool Same(const Trk& a, const Trk& b) { return a.p == b.p && a.th == b.th && a.ph == b.ph; }
int AddDistinct(std::vector<Trk>& v, const Trk& t) {
    for (const auto& x : v) if (Same(x, t)) return (int)v.size();
    v.push_back(t);
    return (int)v.size();
}

struct Stat {
    Long64_t rows = 0, groups = 0, best1 = 0, corbad = 0, grpBestNot1 = 0;
    Long64_t nanEp = 0, nanEm = 0, nanPip = 0, nanPim = 0, nanMee = 0;
    Long64_t fileIdBack = 0;
    double maxd[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};   // recomputed vs stored, cor_ok == 1
    Long64_t nrowHist[6] = {0, 0, 0, 0, 0, 0};          // groups with 1,2,3,4,5,>=6 rows
    Long64_t multiPip = 0, multiPim = 0, multiL1 = 0, multiL2 = 0;
    Long64_t withGam = 0, parallelEE = 0, parallelPP = 0;
    double sumDist[2][2] = {{0, 0}, {0, 0}};            // [convention][lepton/pion], isUsedVertex tracks
    Long64_t nDist[2] = {0, 0};
};

} // namespace prep

void scan_prepare() {
    using namespace prep;
    TStopwatch swAll;
    const std::string prod = scan::kProdDir;
    const char* inName[3] = {"output_pippimepem_exp.root", "output_pippimepep_exp.root",
                             "output_pippimemem_exp.root"};
    const Long64_t expRows[3] = {4783467, 122875, 371524};
    const Long64_t expGrp[3]  = {4058586, 102559, 302701};
    const float nan = std::numeric_limits<float>::quiet_NaN();

    const double mPi = 139.56995, mE = 0.51099895, mP = 938.27231, mEtaMeV = 547.862;
    const double Ebeam = scan::kTbeam + mP;
    const TLorentzVector ini(0, 0, std::sqrt(Ebeam * Ebeam - mP * mP), Ebeam + mP);
    const double yCM = ini.Rapidity();

    TFile* fout = TFile::Open(scan::kCompact, "RECREATE");
    if (!fout || fout->IsZombie()) { printf("cannot create %s\n", scan::kCompact); return; }
    fout->SetCompressionSettings(404);   // LZ4: fast to read

    Stat st[3];
    bool ok = true;

    for (int k = 0; k < 3; ++k) {
        TStopwatch sw;
        TFile* fin = TFile::Open((prod + "/" + inName[k]).c_str(), "READ");
        if (!fin || fin->IsZombie()) { printf("cannot open %s\n", inName[k]); return; }
        TTree* tin = (TTree*)fin->Get("scan_nt");
        if (!tin) { printf("scan_nt missing in %s\n", inName[k]); return; }
        In in;
        tin->SetBranchStatus("*", 0);
#define X(n) tin->SetBranchStatus(#n, 1); tin->SetBranchAddress(#n, &in.n);
        PREP_IN_COLS(X)
#undef X
        tin->SetCacheSize(200 * 1024 * 1024);

        fout->cd();
        TTree* tout = new TTree(scan::SampleName(k),
                                Form("compact scan rows, %s", inName[k]));
        scan::Row row;
        scan::BindRow(tout, row, true);
        tout->SetAutoSave(0);

        std::vector<scan::Row> buf;                    // rows of the current group
        std::vector<Trk> tPip, tPim, tL1, tL2, tLall;
        int nBestInGroup = 0;
        Long64_t grp = -1;
        float keyF = -1, keyR = -1, keyE = -1, lastFile = -1;
        Stat& S = st[k];

        auto flush = [&]() {
            if (buf.empty()) return;
            const int nrow = (int)buf.size();
            for (auto& r : buf) {
                r.nrow = nrow;
                r.n_pip = (int)tPip.size(); r.n_pim = (int)tPim.size();
                r.n_lep1 = (int)tL1.size(); r.n_lep2 = (int)tL2.size();
                r.n_lep = (int)tLall.size();
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
            if (i == 0 || in.file_id != keyF || in.runnumber != keyR || in.event != keyE) {
                flush();
                ++grp;
                keyF = in.file_id; keyR = in.runnumber; keyE = in.event;
                if (in.file_id < lastFile) S.fileIdBack++;
                lastFile = in.file_id;
            }
            scan::Row r;
            // ---- keys ----
            r.file_id = (int)std::lround(in.file_id);
            r.pat_entry = (int)std::lround(in.pat_entry);
            r.run = (int)std::lround(in.runnumber);
            r.grp = (int)grp;
            r.fold = (int)(grp & 1);
            r.isBest = (int)std::lround(in.isBest);
            r.cor_ok = (int)std::lround(in.cor_ok);
            r.nrow = r.n_pip = r.n_pim = r.n_lep1 = r.n_lep2 = r.n_lep = 0;
            if (r.isBest == 1) { nBestInGroup++; S.best1++; }
            if (!r.cor_ok) S.corbad++;
            if (std::isnan(in.ep_p_cor)) S.nanEp++;
            if (std::isnan(in.em_p_cor)) S.nanEm++;
            if (std::isnan(in.pip_p_cor)) S.nanPip++;
            if (std::isnan(in.pim_p_cor)) S.nanPim++;
            if (std::isnan(in.m_ee_cor)) S.nanMee++;

            // ---- stored corrected values ----
            r.m4 = in.m_pippimepem_cor; r.mee = in.m_ee_cor; r.mpipi = in.m_pippim_cor;
            r.oa = in.oa; r.oa_lab = in.oa_pippim_epem_lab_cor;
            r.oa_eta = in.oa_pippim_epem_eta_rest_cor;
            r.mm4 = in.mm_pippimepem_cor; r.mm_pipi = in.mm_pippim_cor;
            r.mm_ee = in.mm_epem_cor; r.pt4 = in.pt4_cor; r.y4 = in.y4_cms_cor;

            // ---- recomputation with the fallback momenta ----
            // corrected momentum where it is a number, reconstructed otherwise
            const double pEp  = std::isnan(in.ep_p_cor)  ? in.ep_p_rec  : in.ep_p_cor;
            const double pEm  = std::isnan(in.em_p_cor)  ? in.em_p_rec  : in.em_p_cor;
            const double pPip = std::isnan(in.pip_p_cor) ? in.pip_p_rec : in.pip_p_cor;
            const double pPim = std::isnan(in.pim_p_cor) ? in.pim_p_rec : in.pim_p_cor;
            const TLorentzVector vEp  = P4(pEp,  in.ep_theta,  in.ep_phi,  mE);
            const TLorentzVector vEm  = P4(pEm,  in.em_theta,  in.em_phi,  mE);
            const TLorentzVector vPip = P4(pPip, in.pip_theta, in.pip_phi, mPi);
            const TLorentzVector vPim = P4(pPim, in.pim_theta, in.pim_phi, mPi);
            const TLorentzVector vEE = vEp + vEm, vPP = vPip + vPim, v4 = vEE + vPP;
            double fb[10];
            fb[0] = v4.M() / 1000.0;
            fb[1] = vEE.M() / 1000.0;
            fb[2] = vPP.M() / 1000.0;
            fb[3] = vEp.Vect().Angle(vEm.Vect()) * kR2D;
            fb[4] = vPP.Vect().Angle(vEE.Vect()) * kR2D;
            {
                const TVector3 ptot = v4.Vect();
                const TLorentzVector hyp(ptot, std::sqrt(ptot.Mag2() + mEtaMeV * mEtaMeV));
                TLorentzVector a = vPP, b = vEE;
                const TVector3 bst = -hyp.BoostVector();
                a.Boost(bst); b.Boost(bst);
                fb[5] = a.Vect().Angle(b.Vect()) * kR2D;
            }
            fb[6] = (ini - v4).M() / 1000.0;
            fb[7] = (ini - vPP).M() / 1000.0;
            fb[8] = v4.Pt();
            fb[9] = v4.Rapidity() - yCM;
            if (r.cor_ok) {
                const double stv[10] = {r.m4, r.mee, r.mpipi, r.oa, r.oa_lab, r.oa_eta,
                                        r.mm4, r.mm_pipi, r.pt4, r.y4};
                for (int j = 0; j < 10; ++j) {
                    const double dd = std::fabs(fb[j] - stv[j]);
                    if (dd > S.maxd[j]) S.maxd[j] = dd;
                }
                r.m4_fb = r.m4; r.mee_fb = r.mee; r.mpipi_fb = r.mpipi; r.oa_fb = r.oa;
                r.oa_lab_fb = r.oa_lab; r.oa_eta_fb = r.oa_eta; r.mm4_fb = r.mm4;
                r.mm_pipi_fb = r.mm_pipi; r.pt4_fb = r.pt4; r.y4_fb = r.y4;
            } else {
                r.m4_fb = (float)fb[0]; r.mee_fb = (float)fb[1]; r.mpipi_fb = (float)fb[2];
                r.oa_fb = (float)fb[3]; r.oa_lab_fb = (float)fb[4]; r.oa_eta_fb = (float)fb[5];
                r.mm4_fb = (float)fb[6]; r.mm_pipi_fb = (float)fb[7];
                r.pt4_fb = (float)fb[8]; r.y4_fb = (float)fb[9];
            }

            // ---- vertex, tracks ----
            r.vx = in.eVertReco_x; r.vy = in.eVertReco_y; r.vz = in.eVertReco_z;
            r.vchi2 = in.eVertReco_chi2;
            r.ep_p = in.ep_p_cor; r.em_p = in.em_p_cor; r.pip_p = in.pip_p_cor; r.pim_p = in.pim_p_cor;
            r.ep_p_rec = in.ep_p_rec; r.em_p_rec = in.em_p_rec;
            r.pip_p_rec = in.pip_p_rec; r.pim_p_rec = in.pim_p_rec;
            r.ep_theta = in.ep_theta; r.em_theta = in.em_theta;
            r.pip_theta = in.pip_theta; r.pim_theta = in.pim_theta;
            r.ep_phi = in.ep_phi; r.em_phi = in.em_phi; r.pip_phi = in.pip_phi; r.pim_phi = in.pim_phi;
            r.ep_rkchi2 = in.ep_rkchi2; r.em_rkchi2 = in.em_rkchi2;
            r.pip_rkchi2 = in.pip_rkchi2; r.pim_rkchi2 = in.pim_rkchi2;
            r.ep_system = (int)std::lround(in.ep_system); r.em_system = (int)std::lround(in.em_system);
            r.pip_system = (int)std::lround(in.pip_system); r.pim_system = (int)std::lround(in.pim_system);
            r.ep_isUsedVertex = (int)std::lround(in.ep_isUsedVertex);
            r.em_isUsedVertex = (int)std::lround(in.em_isUsedVertex);
            r.pip_isUsedVertex = (int)std::lround(in.pip_isUsedVertex);
            r.pim_isUsedVertex = (int)std::lround(in.pim_isUsedVertex);
            r.ep_kIsLepton = (int)std::lround(in.ep_kIsLepton);
            r.em_kIsLepton = (int)std::lround(in.em_kIsLepton);
            r.ep_richmatchquality = in.ep_richmatchquality; r.em_richmatchquality = in.em_richmatchquality;
            r.ep_ringchi2 = in.ep_ringchi2; r.em_ringchi2 = in.em_ringchi2;
            r.ep_oa_lept = in.ep_oa_lept; r.em_oa_lept = in.em_oa_lept;
            r.ep_oa_hadr = in.ep_oa_hadr; r.em_oa_hadr = in.em_oa_hadr;
            r.shared_ring = (in.ep_rich_ind >= 0 && in.ep_rich_ind == in.em_rich_ind) ? 1 : 0;
            r.pip_isring = (int)std::lround(in.pip_isring); r.pim_isring = (int)std::lround(in.pim_isring);
            r.pip_mass2 = in.pip_mass2; r.pim_mass2 = in.pim_mass2;
            r.pip_dedx_mdc = in.pip_dedx_mdc; r.pim_dedx_mdc = in.pim_dedx_mdc;
            r.evtStart = (int)std::lround(in.evtStart);
            r.hyp_chi2 = in.hyp_chi2;
            r.neutr_mult = (int)std::lround(in.neutr_mult);
            r.fwdet_mult = (int)std::lround(in.fwdet_mult);
            r.fw1_p = in.fw1_p; r.fw1_theta = in.fw1_theta;

            // ---- crossing points ----
            const CrossPt cE = Crossing(in.ep_r, in.ep_z, in.ep_theta, in.ep_phi,
                                      in.em_r, in.em_z, in.em_theta, in.em_phi,
                                      in.eVertReco_x, in.eVertReco_y, in.eVertReco_z);
            const CrossPt cP = Crossing(in.pip_r, in.pip_z, in.pip_theta, in.pip_phi,
                                      in.pim_r, in.pim_z, in.pim_theta, in.pim_phi,
                                      in.eVertReco_x, in.eVertReco_y, in.eVertReco_z);
            r.dz_ee = cE.dz; r.r_ee = cE.r; r.rv_ee = cE.rv; r.doca_ee = cE.doca;
            r.dz_pipi = cP.dz; r.r_pipi = cP.r; r.rv_pipi = cP.rv; r.doca_pipi = cP.doca;
            if (cE.r < 0) S.parallelEE++;
            if (cP.r < 0) S.parallelPP++;
            // convention check: tracks used in the vertex fit must pass close to it
            for (int cv = 0; cv < 2; ++cv) {
                if (in.ep_isUsedVertex > 0.5) {
                    S.sumDist[cv][0] += DistToPoint(in.ep_r, in.ep_z, in.ep_theta, in.ep_phi,
                                                    in.eVertReco_x, in.eVertReco_y, in.eVertReco_z, cv);
                    if (cv == 0) S.nDist[0]++;
                }
                if (in.pip_isUsedVertex > 0.5) {
                    S.sumDist[cv][1] += DistToPoint(in.pip_r, in.pip_z, in.pip_theta, in.pip_phi,
                                                    in.eVertReco_x, in.eVertReco_y, in.eVertReco_z, cv);
                    if (cv == 0) S.nDist[1]++;
                }
            }

            // ---- photon with M(e+e-gamma) closest to 0.140 ----
            r.n_gam = (int)std::lround(in.n_gam);
            r.g_idx = -1; r.g_e = r.g_theta = r.g_phi = -1.f;
            r.m_eeg = r.m5 = r.mm5 = r.m5_fb = r.mm5_fb = -1.f;
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
                    r.m5_fb = (float)((v4 + vBest).M() / 1000.0);
                    r.mm5_fb = (float)((ini - v4 - vBest).M() / 1000.0);
                    r.m5 = r.cor_ok ? r.m5_fb : nan;
                    r.mm5 = r.cor_ok ? r.mm5_fb : nan;
                }
            }

            // ---- derived ----
            r.T_pi0 = (float)scan::TPi0(r.mpipi);        // NaN when cor_ok == 0
            r.m4max = std::isnan(r.mpipi) ? nan : (float)scan::M4max(r.mpipi);

            // ---- distinct tracks of the group ----
            const Trk a = {in.pip_p_rec, in.pip_theta, in.pip_phi};
            const Trk b = {in.pim_p_rec, in.pim_theta, in.pim_phi};
            const Trk c = {in.ep_p_rec, in.ep_theta, in.ep_phi};
            const Trk d = {in.em_p_rec, in.em_theta, in.em_phi};
            AddDistinct(tPip, a); AddDistinct(tPim, b);
            AddDistinct(tL1, c); AddDistinct(tL2, d);
            AddDistinct(tLall, c); AddDistinct(tLall, d);

            buf.push_back(r);
            S.rows++;
        }
        flush();
        S.groups = grp + 1;

        fout->cd();
        tout->Write("", TObject::kOverwrite);
        fin->Close();
        delete fin;
        printf("[%s] %lld rows, %lld groups, %.1f s\n", scan::SampleName(k), S.rows, S.groups,
               sw.RealTime());
        if (S.rows != expRows[k] || S.groups != expGrp[k] || S.best1 != expGrp[k]) ok = false;
    }
    fout->Close();
    delete fout;

    // ---- read-back check: compact file against the source, row by row ----
    Long64_t rbRows[3] = {0, 0, 0}, rbGroups[3] = {0, 0, 0}, rbBest[3] = {0, 0, 0};
    Long64_t rbM4bad[3] = {0, 0, 0}, rbCorOk[3] = {0, 0, 0}, rbKeyBad[3] = {0, 0, 0};
    {
        TFile* fc = TFile::Open(scan::kCompact, "READ");
        for (int k = 0; k < 3 && fc && !fc->IsZombie(); ++k) {
            TTree* tc = (TTree*)fc->Get(scan::SampleName(k));
            TFile* fin = TFile::Open((prod + "/" + inName[k]).c_str(), "READ");
            TTree* tin = fin ? (TTree*)fin->Get("scan_nt") : nullptr;
            if (!tc || !tin) { ok = false; continue; }
            Int_t cGrp = 0, cBest = 0, cCor = 0, cPat = 0, cFile = 0;
            Float_t cM4 = 0, sM4 = 0, sPat = 0, sFile = 0;
            tc->SetBranchStatus("*", 0);
            for (const char* b : {"grp", "isBest", "cor_ok", "pat_entry", "file_id", "m4"})
                tc->SetBranchStatus(b, 1);
            tc->SetBranchAddress("grp", &cGrp); tc->SetBranchAddress("isBest", &cBest);
            tc->SetBranchAddress("cor_ok", &cCor); tc->SetBranchAddress("pat_entry", &cPat);
            tc->SetBranchAddress("file_id", &cFile); tc->SetBranchAddress("m4", &cM4);
            tin->SetBranchStatus("*", 0);
            for (const char* b : {"m_pippimepem_cor", "pat_entry", "file_id"}) tin->SetBranchStatus(b, 1);
            tin->SetBranchAddress("m_pippimepem_cor", &sM4);
            tin->SetBranchAddress("pat_entry", &sPat); tin->SetBranchAddress("file_id", &sFile);
            const Long64_t n = tc->GetEntries();
            rbRows[k] = n;
            if (n != tin->GetEntries()) ok = false;
            Int_t lastGrp = -1;
            for (Long64_t i = 0; i < n; ++i) {
                tc->GetEntry(i); tin->GetEntry(i);
                if (cGrp != lastGrp) { rbGroups[k]++; lastGrp = cGrp; }
                if (cBest == 1) rbBest[k]++;
                if (cPat != (int)std::lround(sPat) || cFile != (int)std::lround(sFile)) rbKeyBad[k]++;
                if (cCor) { rbCorOk[k]++; if (!(cM4 == sM4)) rbM4bad[k]++; }
            }
            if (rbRows[k] != expRows[k] || rbGroups[k] != expGrp[k] || rbBest[k] != expGrp[k] ||
                rbM4bad[k] != 0 || rbKeyBad[k] != 0) ok = false;
            fin->Close();
            delete fin;
        }
        if (fc) { fc->Close(); delete fc; }
    }

    // ---- verification report ----
    const std::string logName = std::string(scan::kScanDir) + "/scan_prepare.log";
    FILE* lf = fopen(logName.c_str(), "w");
    const char* vn[10] = {"m4", "mee", "mpipi", "oa", "oa_lab", "oa_eta", "mm4", "mm_pipi", "pt4", "y4"};
    for (int pass = 0; pass < 2; ++pass) {
        FILE* o = pass == 0 ? stdout : lf;
        if (!o) continue;
        fprintf(o, "scan_prepare: verification of %s\n", scan::kCompact);
        for (int k = 0; k < 3; ++k) {
            const Stat& S = st[k];
            fprintf(o, "\n[%s] %s\n", scan::SampleName(k), inName[k]);
            fprintf(o, "  rows %lld (expected %lld)  groups %lld  isBest==1 rows %lld (expected %lld)\n",
                    S.rows, expRows[k], S.groups, S.best1, expGrp[k]);
            fprintf(o, "  groups without exactly one isBest==1 row: %lld ; file_id decreasing: %lld times\n",
                    S.grpBestNot1, S.fileIdBack);
            fprintf(o, "  cor_ok==0 rows %lld ; NaN p_cor: e1 %lld, e2 %lld, pi+ %lld, pi- %lld ; NaN m_ee_cor %lld\n",
                    S.corbad, S.nanEp, S.nanEm, S.nanPip, S.nanPim, S.nanMee);
            fprintf(o, "  read-back of the compact tree: rows %lld, groups %lld, isBest==1 %lld, cor_ok==1 %lld\n",
                    rbRows[k], rbGroups[k], rbBest[k], rbCorOk[k]);
            fprintf(o, "  read-back: m4 != m_pippimepem_cor for cor_ok==1 in %lld rows ; (file_id, pat_entry) differ in %lld rows\n",
                    rbM4bad[k], rbKeyBad[k]);
            fprintf(o, "  max |recomputed - stored| for cor_ok==1:");
            for (int j = 0; j < 10; ++j) fprintf(o, " %s %.2g", vn[j], S.maxd[j]);
            fprintf(o, "\n  groups with 1,2,3,4,5,>=6 rows: %lld %lld %lld %lld %lld %lld\n",
                    S.nrowHist[0], S.nrowHist[1], S.nrowHist[2], S.nrowHist[3], S.nrowHist[4], S.nrowHist[5]);
            fprintf(o, "  groups with more than one distinct pi+ %lld, pi- %lld, lepton1 %lld, lepton2 %lld\n",
                    S.multiPip, S.multiPim, S.multiL1, S.multiL2);
            fprintf(o, "  rows with a photon %lld ; parallel lines: leptons %lld, pions %lld\n",
                    S.withGam, S.parallelEE, S.parallelPP);
            fprintf(o, "  mean distance track line - event vertex [mm], tracks used in the vertex fit:\n");
            fprintf(o, "    HYDRA convention (base = r*(-sin phi, cos phi)): lepton1 %.3f, pi+ %.3f\n",
                    S.nDist[0] ? S.sumDist[0][0] / S.nDist[0] : 0., S.nDist[1] ? S.sumDist[0][1] / S.nDist[1] : 0.);
            fprintf(o, "    opposite sign of r:                              lepton1 %.3f, pi+ %.3f\n",
                    S.nDist[0] ? S.sumDist[1][0] / S.nDist[0] : 0., S.nDist[1] ? S.sumDist[1][1] / S.nDist[1] : 0.);
        }
        fprintf(o, "\nresult: %s\n", ok ? "row and group counts as expected" : "MISMATCH in row or group counts");
    }
    if (lf) fclose(lf);
    printf("scan_prepare: total %.1f s, log %s\n", swAll.RealTime(), logName.c_str());
}
