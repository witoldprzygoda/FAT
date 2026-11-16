// ============================================================================
// sim_eta_physics.h  --  physics reweighting of the simulated eta decays
// ============================================================================
// HADES p+p, T = 4.53 GeV, pi+ pi- e+ e-. Weights (no new production) that bring
// the eta decays of the SMASH simulation to physics models:
//
//   direct eta -> pi+pi-e+e- (four-particle class 5, eta set; generator: PLUTO
//     "eta:e+,e-,pi+,pi-" added with makeStaticData()->AddDecay, eta_physics
//     plugin disabled, NOT QED in M(ee)):
//       w5 = T(q, M[, angles]) / G(q, M[, angles]),  q = M(ee)_true, M = M(pipi)_true
//   class 7 of the eta set and every other eta-set row with the e+e- pair of a
//     forced decay (lepton-pair class 10): M(ee)-only weight T_ee(q) / G_ee(q)
//     (also in V3/V4: their e+ angles are not in the angle file, UNVERIFIED approximation);
//     eta-set rows with the pion pair of a forced decay (pion class 9) but not
//     class 5: M(pipi)-only weight T_pp(M) / G_pp(M)   (marginals of T and G)
//   X17 (class 6): unweighted
//   bump, eta -> pi+pi-pi0 (classes 2, 3; HGeant, phase space ASSUMED):
//       wB = f(X, Y) / <f>_PS with the X^2 terms replaced by their phase-space
//       average at fixed Y (Y from M(pipi)_true; X unknown: the photon is lost)
//   peak conversion eta -> pi+pi-gamma (class 4; HGeant, phase space ASSUMED):
//       wC = |F_V|^2 p*^2 k^2 / <|F_V|^2 p*^2 k^2>_PS
//   other rows whose two pions come from one eta (pion class 2, any set): the
//       HGeant mixture of pi+pi-pi0 and pi+pi-gamma at the row's M(pipi)_true,
//       w = p3(M) wB(M) + (1 - p3(M)) wC(M), p3 = B3 rho3 / (B3 rho3 + Bg rhog)
//       with the HGeant branching ratios B3 = 0.2114, Bg = 0.0478 (sim_truth.h)
//       and the phase-space M(pipi) densities rho of the two decays (pion class 2
//       does not tell the two modes apart; above m_eta - m_pi0 only pi+pi-gamma)
//
// TARGET MODEL T of eta -> pi+pi-e+e- (density in q = M(ee) and M = M(pipi)):
//   T(q, M) = KW(q) x R(q, M)
//   KW(q)   = (2 alpha / (3 pi q)) (1 + 2 m_e^2/q^2) sqrt(1 - 4 m_e^2/q^2)
//             [Kroll-Wada conversion of a transverse virtual photon into e+e-,
//              N. M. Kroll, W. Wada, Phys. Rev. 98 (1955) 1355; per dq]
//   R(q, M) = p*^3(M) lambda^{3/2}(m_eta^2, M^2, q^2) |F_V(M^2)|^2 |F_VMD(q^2)|^2
//             [rate dGamma(eta -> pi+pi- gamma*)/dM of the anomalous P-wave amplitude
//              A ~ eps_{mu nu a b} eps^mu p+^nu p-^a k^b: |A|^2 ~ p*^2 lambda sin^2(theta_pi),
//              phase space ~ (lambda^{1/2} / m_eta^2)(p* / M) dM^2, so dGamma/dM ~ p*^3 lambda^{3/2};
//              p* = pion momentum in the pi+pi- rest frame = M beta_pi / 2, i.e. the factor is
//              (beta_pi^3 M^3 / 8) per dM = (s beta_pi^3 / 8) per ds; for q = 0, lambda^{1/2} =
//              m_eta^2 - M^2 = 2 m_eta k: the real-photon rate p*^3 k^3 |F_V|^2 (e.g. Stollenwerk et
//              al., PLB 707 (2012) 184, dGamma/ds ~ s sigma^3 (m_eta^2 - s)^3 |A|^2)]
//   F_V     pion vector form factor, default the resummed effective-field-theory form of
//           F. Guerrero, A. Pich, PLB 412 (1997) 382, hep-ph/9707347 (pi and K loops,
//           f_pi = 92.4 MeV; Gamma_rho(s) scaled to Gamma_rho(m_rho^2) = 0.1491 GeV);
//           variant: Gounaris-Sakurai, PRL 21 (1968) 244; m_rho = 0.7753, Gamma_rho = 0.1491 GeV
//   F_VMD   rho-pole transition form factor m_V^2 / (m_V^2 - q^2 - i m_V Gamma(q^2)),
//           m_V = m_rho; Gamma(q^2) = 0 below 4 m_pi^2, i.e. in the whole range of q
//           (q < m_eta - 2 m_pi = 0.269 GeV): |F_VMD|^2 = (m_V^2 / (m_V^2 - q^2))^2
//   angles  (variant V3) (1 + cos^2 th_e + (4 m_e^2/q^2) sin^2 th_e) for the e+ in the e+e-
//           rest frame w.r.t. the e+e- direction in the eta frame (transverse gamma*), and
//           sin^2 th_pi for the pi+ in the pi+pi- rest frame w.r.t. the pi+pi- direction
//           (P wave); the generator angular densities are measured (sim_rw_build.C:
//           e+ flat, pi+ sin^2 already: the pion factor is then 1). These factors follow
//           from the same amplitude as p*^3 lambda^{3/2}: in the gamma* rest frame the
//           hadronic current is along p+ x p- (normal to the pion plane, transverse to the
//           gamma* direction), so the full angular density is
//             sin^2 th_pi [1 - beta_e^2 sin^2 th_e sin^2 phi],  beta_e^2 = 1 - 4 m_e^2/q^2,
//           phi = angle between the e+e- and pi+pi- decay planes; its average over phi is the
//           e+ factor above (V3).
//   plane term (variant V4, check only) R_phi = 1 + beta_e^2 sin^2 th_e cos 2phi / (2 - beta_e^2 sin^2 th_e)
//           = full density / V3 density; generator phi ASSUMED flat (PLUTO; the generator
//           file has no phi), so R_phi keeps the generator-level mean; checked on the accepted
//           rows in sim_rw_build.C (tables/rw_weights_planes.txt): below 0.5 %
//   variants V0 = KW x p*^3 lambda^{3/2} (no form factors, no angles); V1 = V0 x |F_V|^2;
//           V2 = V1 x |F_VMD|^2; V3 = V2 x angular factors (MAIN since the review of
//           2026-11-12: V2 leaves the generator's flat e+ angle in place, which is not the
//           physics of the same amplitude); V4 = V3 x R_phi (check of the plane term)
//   T and G are each normalised to unit integral over the full phase space
//   (q >= 2 m_e, M >= 2 m_pi, q + M <= m_eta), then the generator-level mean of
//   w5 is set to 1 exactly by the factor 'renorm' measured on the generator sample
//   (it differs from 1 only through the finite description of G), so the number
//   of generated decays (BR 2.68e-4) is unchanged.
//   Threshold bin: the generator file has a limited momentum precision (4.9 % of
//   its M(ee) values lie below 2 m_e); below q = 1.553 MeV (first bin of G) both T
//   and G are used as bin averages in q.
//
// GENERATOR DENSITY G (sim_rw_build.C, file kRwGenFile): factorised
//   G(q, M) = g(q) h(M) lambda^{1/2}(m_eta^2, M^2, q^2), m_eta = 0.548 GeV (PLUTO),
//   g, h fitted by maximum likelihood as bin contents (q: log bins, M: 2 MeV bins) to
//   1 181 801 generated decays (the lambda^{1/2} factorisation is tested there); h is used
//   piecewise constant, g log-log interpolated between bin centres (constant in the
//   threshold bin), see Gq().
//
// BUMP: f(X,Y) = 1 + aY + bY^2 + dX^2 + fY^3 + gX^2Y, KLOE-2 2016 (A. Anastasi et al.,
//   JHEP 05 (2016) 019, arXiv:1601.06985): a = -1.095, b = 0.145, d = 0.081,
//   f = 0.141, g = -0.044; X = sqrt(3)(T+ - T-)/Q, Y = 3 T0/Q - 1, Q = T+ + T- + T0,
//   T0 = (m_eta^2 + m_pi0^2 - M^2)/(2 m_eta) - m_pi0. At fixed Y the phase space is
//   uniform in X between +-Xmax(Y), Xmax = 2 sqrt(3) p* q0* / (m_eta Q) (p*, q0* =
//   pi+- and pi0 momenta in the pi+pi- rest frame), so <X^2>_Y = Xmax^2 / 3.
//   Variant 'slides' (K. Proscinski, X17 analysis): f = 1 + aY, a = -1.095.
//   <f>_PS from TGenPhaseSpace (sim_rw_build.C), stored in kRwGenFile.
//
// USE: simphys::Init() reads kRwGenFile; simphys::Install() connects the weight to
// SimEngine.h (SimOpt::phys, default 0 = OFF: every existing result unchanged).
//   SimOpt o; o.phys = simphys::kPhysMain;   // V3 + KLOE bump + conversion + background rules
// (kPhysMain = 108 since the review of 2026-11-12; the earlier main V2 is kPhysV2 = 107)
// The variant code: bits 0-2 direct (0 off, 1 V0, 2 V1, 3 V2, 4 V3, 5 V4), bits 3-4 bump
// (0 off, 1 KLOE, 2 slides 1 + aY), bit 5 conversion, bit 6 background rules (class 7 /
// lepton class 10 / pion class 9 of the eta set, pion class 2), bit 7 F_V Gounaris-Sakurai.
// V3 and V4 need the true angles: simphys::LoadAngles() (file kRwAngFile, built by
// sim_rw_build.C from the true momenta in scan_nt; key = event group + true momenta);
// RowFactor loads them itself (once) if a V3/V4 code is used before LoadAngles().
// All arithmetic in double.
// ============================================================================
#ifndef SIM_ETA_PHYSICS_H
#define SIM_ETA_PHYSICS_H

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include "SimEngine.h"

namespace simphys {

using simeng::Row;

// ---------------------------------------------------------------- files
const char* const kRwGenFile = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/rw_generator.root";
const char* const kRwAngFile = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/rw_angles_eta.root";
const char* const kGenText =
    "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/generator/eta_pipiee_gen_kin_1000files.txt.gz";

// ---------------------------------------------------------------- constants (GeV)
const double kAlpha = 1. / 137.035999084;
const double kMe = 0.00051099895, kMpi = 0.13957039, kMpi0 = 0.1349768, kMK = 0.495644;   // PDG; m_K isospin average
const double kFpi = 0.0924;                         // f_pi, Guerrero-Pich
const double kMrho = 0.7753, kGrho = 0.1491;        // rho mass and width (as in the X17 slides)
const double kMetaGen = 0.548;                      // PLUTO eta mass of the eta set (generator file, column 5)
const double kMetaHG = 0.547862;                    // eta mass for the HGeant decays (PDG; HGeant's own value UNVERIFIED)
const double kBRppee = 2.68e-4, kBRppg = 4.22e-2;   // PDG eta -> pi+pi-e+e-, pi+pi-gamma
// KLOE-2 2016 eta -> pi+pi-pi0 Dalitz-plot parameters (JHEP 05 (2016) 019)
const double kDalA = -1.095, kDalB = 0.145, kDalD = 0.081, kDalF = 0.141, kDalG = -0.044;
// HGeant branching ratios (src/sim_truth.h) for the pion-class-2 mixture
const double kHGppp0 = simtruth::kHG_etaPPP0, kHGppg = simtruth::kHG_etaPPG;

// ---------------------------------------------------------------- variant code
enum { kDirOff = 0, kV0 = 1, kV1 = 2, kV2 = 3, kV3 = 4, kV4 = 5 };
enum { kBumpOff = 0, kBumpKLOE = 1, kBumpSlides = 2 };
inline int Code(int dir, int bump, int conv, int bkg, int gs = 0) {
    return (dir & 7) | ((bump & 3) << 3) | ((conv & 1) << 5) | ((bkg & 1) << 6) | ((gs & 1) << 7);
}
inline int CodeDir(int c) { return c & 7; }
inline int CodeBump(int c) { return (c >> 3) & 3; }
inline int CodeConv(int c) { return (c >> 5) & 1; }
inline int CodeBkg(int c) { return (c >> 6) & 1; }
inline int CodeGS(int c) { return (c >> 7) & 1; }
const int kPhysMain = 4 | (1 << 3) | (1 << 5) | (1 << 6);   // 108: V3 + KLOE bump + conversion + background rules (MAIN)
const int kPhysV2 = 3 | (1 << 3) | (1 << 5) | (1 << 6);     // 107: V2 (no decay angles) + KLOE + conversion + bkg rules (main before the review)
inline std::string PhysLabel(int c) {
    if (!c) return "no physics weights";
    static const char* dn[6] = {"direct off", "direct V0", "direct V1", "direct V2", "direct V3", "direct V4 (V3 x plane term)"};
    static const char* bn[3] = {"bump off", "bump KLOE", "bump 1+aY"};
    std::string s = std::string(dn[CodeDir(c) <= 5 ? CodeDir(c) : 0]) + ", " + bn[CodeBump(c) <= 2 ? CodeBump(c) : 0] +
                    (CodeConv(c) ? ", conv F_V" : ", conv off") + (CodeBkg(c) ? ", bkg rules" : ", bkg off");
    if (CodeGS(c)) s += ", F_V GS";
    return s;
}

// ---------------------------------------------------------------- kinematics
inline double Lam(double a, double b, double c) { return a * a + b * b + c * c - 2 * a * b - 2 * a * c - 2 * b * c; }
inline double PStar(double M, double m) { const double x = 0.25 * M * M - m * m; return x > 0 ? std::sqrt(x) : 0.; }
// momentum of the third body in the eta frame (eta -> (M) + m3)
inline double K3(double meta, double M, double m3) {
    const double l = Lam(meta * meta, M * M, m3 * m3);
    return l > 0 ? std::sqrt(l) / (2 * meta) : 0.;
}

// ---------------------------------------------------------------- QED
// Kroll-Wada density per dq [1/GeV], relative to the real-photon rate
inline double KW(double q) {
    if (q <= 2 * kMe) return 0.;
    const double r = kMe * kMe / (q * q);
    return (2 * kAlpha / (3 * M_PI * q)) * (1 + 2 * r) * std::sqrt(1 - 4 * r);
}
// normalised e+ angular density per dcos(theta_e) (transverse virtual photon)
inline double AngE(double c, double q) {
    const double r = (q > 2 * kMe) ? kMe * kMe / (q * q) : 0.25;
    return 0.375 * (1 + c * c + 4 * r * (1 - c * c)) / (1 + 2 * r);
}
// normalised pi+ angular density per dcos(theta_pi) (P wave)
inline double AngPi(double c) { return 0.75 * (1 - c * c); }

// ---------------------------------------------------------------- form factors
// Guerrero-Pich loop function: Re A(m^2/s, m^2/mu^2) = ln(m^2/mu^2) + 8 m^2/s - 5/3 + Re[sigma^3 ln((sigma+1)/(sigma-1))]
inline double GPReA(double s, double m) {
    if (s <= 0) return std::log(m * m / (kMrho * kMrho)) + 1.;   // s -> 0 limit (finite)
    const double x = 4 * m * m / s;
    double t;
    if (x < 1) { const double sg = std::sqrt(1 - x); t = sg * sg * sg * std::log((1 + sg) / (1 - sg)); }
    else { const double a = std::sqrt(x - 1); t = -2 * a * a * a * std::atan(1. / a); }   // below threshold
    return std::log(m * m / (kMrho * kMrho)) + x * 2. - 5. / 3. + t;                    // 8 m^2/s = 2 x
}
inline double SigmaCube(double s, double m) { const double x = 1 - 4 * m * m / s; return x > 0 ? x * std::sqrt(x) : 0.; }
// F_V, Guerrero-Pich (hep-ph/9707347); width scaled to Gamma_rho at m_rho^2
inline std::complex<double> FV_GP(double s) {
    const double mr2 = kMrho * kMrho;
    const double sc = SigmaCube(mr2, kMpi);
    const double gam = kGrho * (s / mr2) * (SigmaCube(s, kMpi) + 0.5 * SigmaCube(s, kMK)) / sc;
    const std::complex<double> bw = mr2 / std::complex<double>(mr2 - s, -kMrho * gam);
    const double ex = -(s / (96 * M_PI * M_PI * kFpi * kFpi)) * (GPReA(s, kMpi) + 0.5 * GPReA(s, kMK));
    return bw * std::exp(ex);
}
// F_V, Gounaris-Sakurai (PRL 21 (1968) 244)
inline double GSk(double s) { return 0.5 * std::sqrt(s) * std::sqrt(std::max(0., 1 - 4 * kMpi * kMpi / s)); }
inline double GSh(double s) {
    const double k = GSk(s), rs = std::sqrt(s);
    return (2. / M_PI) * (k / rs) * std::log((rs + 2 * k) / (2 * kMpi));
}
inline std::complex<double> FV_GS(double s) {
    const double mr2 = kMrho * kMrho, km = GSk(mr2);
    const double eps = 1e-6;
    const double hp = (GSh(mr2 * (1 + eps)) - GSh(mr2 * (1 - eps))) / (2 * eps * mr2);
    const double d = (3. / M_PI) * kMpi * kMpi / (km * km) * std::log((kMrho + 2 * km) / (2 * kMpi)) + kMrho / (2 * M_PI * km) -
                     kMpi * kMpi * kMrho / (M_PI * km * km * km);
    double f = 0., gam = 0.;
    if (s > 4 * kMpi * kMpi) {
        const double k = GSk(s);
        f = kGrho * mr2 / (km * km * km) * (k * k * (GSh(s) - GSh(mr2)) + (mr2 - s) * km * km * hp);
        gam = kGrho * (kMrho / std::sqrt(s)) * std::pow(k / km, 3);
    } else {
        // below threshold (not used: s >= 4 m_pi^2 here); analytic continuation omitted
        f = kGrho * mr2 / (km * km * km) * ((mr2 - s) * km * km * hp - km * km * 0.);
    }
    return std::complex<double>(mr2 + d * kMrho * kGrho, 0.) / std::complex<double>(mr2 - s + f, -kMrho * gam);
}
inline double FV2(double s, int gs) { return std::norm(gs ? FV_GS(s) : FV_GP(s)); }
inline double FVMD2(double q2) { const double mv2 = kMrho * kMrho; return (mv2 / (mv2 - q2)) * (mv2 / (mv2 - q2)); }

// ---------------------------------------------------------------- bump: eta -> pi+pi-pi0 Dalitz plot
inline double DalQ(double meta) { return meta - 2 * kMpi - kMpi0; }
inline double DalY(double M, double meta = kMetaHG) {
    const double t0 = (meta * meta + kMpi0 * kMpi0 - M * M) / (2 * meta) - kMpi0;
    return 3 * t0 / DalQ(meta) - 1;
}
inline double DalXmax(double M, double meta = kMetaHG) {
    const double ps = PStar(M, kMpi);
    const double l = Lam(meta * meta, M * M, kMpi0 * kMpi0);
    const double q0 = (l > 0 && M > 0) ? std::sqrt(l) / (2 * M) : 0.;
    return 2 * std::sqrt(3.) * ps * q0 / (meta * DalQ(meta));
}
inline double DalF(double X, double Y, int bv) {
    if (bv == kBumpSlides) return 1 + kDalA * Y;
    return 1 + kDalA * Y + kDalB * Y * Y + kDalD * X * X + kDalF * Y * Y * Y + kDalG * X * X * Y;
}
// f at fixed Y with X^2 -> <X^2>_Y = Xmax^2 / 3 (phase space uniform in X at fixed Y)
inline double DalFbar(double M, int bv, double meta = kMetaHG) {
    const double Mc = std::min(std::max(M, 2 * kMpi), meta - kMpi0);
    const double Y = DalY(Mc, meta);
    if (bv == kBumpSlides) return 1 + kDalA * Y;
    const double x2 = DalXmax(Mc, meta) * DalXmax(Mc, meta) / 3.;
    return 1 + kDalA * Y + kDalB * Y * Y + kDalF * Y * Y * Y + (kDalD + kDalG * Y) * x2;
}
// HGeant phase-space M(pipi) densities (unnormalised): eta -> pi+pi- X, dPhi3/dM ~ p* k_X
inline double Rho3(double M, double meta = kMetaHG) { return (M > 2 * kMpi && M < meta - kMpi0) ? PStar(M, kMpi) * K3(meta, M, kMpi0) : 0.; }
inline double RhoG(double M, double meta = kMetaHG) { return (M > 2 * kMpi && M < meta) ? PStar(M, kMpi) * K3(meta, M, 0.) : 0.; }

// ---------------------------------------------------------------- the model (singleton)
struct Model {
    bool ok = false;
    // generator density G(q, M) = g(q) h(M) lambda^{1/2}
    std::vector<double> qE, g, mE, h;
    // fine tables
    int nQ = 0, nM = 0;
    double lq0 = 0, lq1 = 0, m0 = 0, m1 = 0;
    std::vector<double> Hlam;        // H(q) = int h lambda^{1/2} dM on the q grid
    std::vector<double> Glam;        // int g lambda^{1/2} dq on the M grid
    std::vector<double> Iee[3][2];   // int R_v dM (mass part, no KW) on the q grid, v = 0,1,2, fv = GP/GS
    std::vector<double> Tpp[3][2];   // int T_v dq on the M grid (unnormalised)
    std::vector<double> fvTab[2];    // |F_V|^2 on a fine M grid
    double fvM0 = 0, fvM1 = 0; int fvN = 0;
    double ZT[3][2] = {{1, 1}, {1, 1}, {1, 1}};
    double KWbar0 = 0;               // bin average of KW over the threshold bin of g
    // renormalisation measured on the generator sample (1 = none): W5 [v=0..3][fv], Wee [v=0..2][fv], Wpp [v][fv]
    double ren5[4][2] = {{1, 1}, {1, 1}, {1, 1}, {1, 1}}, renEE[3][2] = {{1, 1}, {1, 1}, {1, 1}}, renPP[3][2] = {{1, 1}, {1, 1}, {1, 1}};
    // generator angular densities: a(c) = (1 + b c^2) / (2 + 2 b / 3); bPi = -1 means sin^2 (pion factor 1)
    double bE = 0., bPi = -1.;
    // bump and conversion normalisations
    double normB[3] = {1, 1, 1};     // <f>_PS (index kBumpKLOE, kBumpSlides), TGenPhaseSpace
    double normC[2] = {1, 1};        // <|F_V|^2 p*^2 k^2>_PS of eta -> pi+pi-gamma (GP, GS)
    double rho3N = 1, rhoGN = 1;     // normalisations of Rho3, RhoG
};
inline Model& Mdl() { static Model m; return m; }

// helpers on the model
inline int BinOf(const std::vector<double>& e, double x) {   // clamp to [0, n-1]
    const int n = (int)e.size() - 1;
    if (x < e[1]) return 0;
    if (x >= e[n - 1]) return n - 1;
    int lo = 1, hi = n - 1;
    while (hi - lo > 1) { const int mid = (lo + hi) / 2; if (x < e[mid]) hi = mid; else lo = mid; }
    return lo;
}
inline double Interp(const std::vector<double>& t, double x0, double x1, double x) {
    const int n = (int)t.size();
    double u = (x - x0) / (x1 - x0) * (n - 1);
    if (u <= 0) return t[0];
    if (u >= n - 1) return t[n - 1];
    const int i = (int)u; const double f = u - i;
    return t[i] * (1 - f) + t[i + 1] * f;
}
inline double FV2t(double x, int gs) {   // tabulated |F_V(x^2)|^2, x = M(pipi)
    const Model& m = Mdl();
    if (m.fvN == 0 || x < m.fvM0 || x > m.fvM1) return FV2(x * x, gs);
    return Interp(m.fvTab[gs], m.fvM0, m.fvM1, x);
}
inline double QMax() { return kMetaGen - 2 * kMpi; }
// g(q) of the generator density: constant in the threshold bin [2 m_e, qE[1]); above, log-log interpolation between the
// geometric bin centres (bin contents of the maximum-likelihood fit), extrapolated with the slope of the first two bins
// below the first centre and constant above the last centre (removes the saw-tooth of a piecewise-constant g, whose
// in-bin slope (about q^-2) would otherwise survive in the reweighted sample)
inline double Gq(double q) {
    const Model& m = Mdl();
    const int n = (int)m.g.size();
    if (q < m.qE[1]) return m.g[0];
    const double lq = std::log(q);
    auto lc = [&](int a) { return 0.5 * (std::log(m.qE[a]) + std::log(m.qE[a + 1])); };
    if (lq >= lc(n - 1)) return m.g[n - 1];
    int a = 1;
    if (lq < lc(1)) a = 1;
    else { a = BinOf(m.qE, q); if (lq < lc(a)) --a; }
    if (a >= n - 1) a = n - 2;
    if (!(m.g[a] > 0) || !(m.g[a + 1] > 0)) return m.g[BinOf(m.qE, q)];
    const double t = (lq - lc(a)) / (lc(a + 1) - lc(a));
    return std::exp(std::log(m.g[a]) * (1 - t) + std::log(m.g[a + 1]) * t);
}
inline double KWeff(double q) {          // KW, as a bin average in the threshold bin of g
    const Model& m = Mdl();
    if (!m.qE.empty() && q < m.qE[1]) return m.KWbar0;
    return KW(q);
}
// mass part of T without KW (p*^3 lambda^{3/2} FF), and lambda^{1/2}
inline double Rmass(double q, double Mpp, int v, int gs) {
    const double l = Lam(kMetaGen * kMetaGen, Mpp * Mpp, q * q);
    if (l <= 0 || Mpp <= 2 * kMpi) return 0.;
    const double ps = PStar(Mpp, kMpi);
    double r = ps * ps * ps * l * std::sqrt(l);
    if (v >= 1) r *= FV2t(Mpp, gs);
    if (v >= 2) r *= FVMD2(q * q);
    return r;
}
// generator density G (normalised)
inline double Gdens(double q, double Mpp) {
    const Model& m = Mdl();
    const double l = Lam(kMetaGen * kMetaGen, Mpp * Mpp, q * q);
    if (l <= 0) return 0.;
    return Gq(q) * m.h[BinOf(m.mE, Mpp)] * std::sqrt(l);
}
// target density T_v (normalised; v = 0, 1, 2)
inline double Tdens(double q, double Mpp, int v, int gs) { return KWeff(q) * Rmass(q, Mpp, v, gs) / Mdl().ZT[v][gs]; }

// plane term R_phi = full angular density / V3 density (see the header comment); c2p = cos(2 phi)
inline double RPhi(double ce, double c2p, double q) {
    const double b2 = (q > 2 * kMe) ? 1 - 4 * kMe * kMe / (q * q) : 0.;
    const double x = b2 * (1 - ce * ce);
    return 1 + x * c2p / (2 - x);
}
// direct-decay weight (class 5); dv = kV0..kV4; ce, cp = true cos(theta*) of e+ and pi+, c2p = cos(2 phi) (V3, V4)
inline double W5(double q, double Mpp, int dv, int gs, double ce = 0, double cp = 0, bool haveAng = false, double c2p = 0) {
    const Model& m = Mdl();
    if (dv == kDirOff) return 1.;
    const int v = dv >= kV3 ? 2 : dv - 1;
    q = std::max(q, 2 * kMe * (1 + 1e-9));
    const double l = Lam(kMetaGen * kMetaGen, Mpp * Mpp, q * q);
    if (l <= 0 || Mpp <= 2 * kMpi) return 0.;
    const double ps = PStar(Mpp, kMpi);
    double r = ps * ps * ps * l;                     // lambda^{3/2} / lambda^{1/2}
    if (v >= 1) r *= FV2t(Mpp, gs);
    if (v >= 2) r *= FVMD2(q * q);
    const double gq = Gq(q), hm = m.h[BinOf(m.mE, Mpp)];
    if (!(gq > 0) || !(hm > 0)) return 0.;
    double w = KWeff(q) * r / m.ZT[v][gs] / (gq * hm);
    if (dv >= kV3) {
        w *= m.ren5[3][gs] / m.ren5[2][gs];
        if (haveAng) {
            const double ae = (1 + m.bE * ce * ce) / (2 + 2 * m.bE / 3);
            w *= AngE(ce, q) / ae;
            if (m.bPi != -1.) {
                const double ap = (1 + m.bPi * cp * cp) / (2 + 2 * m.bPi / 3);
                w *= AngPi(cp) / ap;
            }
            if (dv == kV4) w *= RPhi(ce, c2p, q);
        }
    }
    return w * m.ren5[v][gs];
}
// M(ee)-only weight (lepton class 10 rows other than class 5): T_ee / G_ee
inline double Wee(double q, int dv, int gs) {
    const Model& m = Mdl();
    if (dv == kDirOff) return 1.;
    const int v = dv >= kV3 ? 2 : dv - 1;
    q = std::max(q, 2 * kMe * (1 + 1e-9));
    if (q >= QMax()) return 0.;
    const double lq = std::log(q);
    const double te = KWeff(q) * Interp(m.Iee[v][gs], m.lq0, m.lq1, lq) / m.ZT[v][gs];
    const double ge = Gq(q) * Interp(m.Hlam, m.lq0, m.lq1, lq);
    return ge > 0 ? m.renEE[v][gs] * te / ge : 0.;
}
// M(pipi)-only weight (pion class 9 rows other than class 5): T_pp / G_pp
inline double Wpp(double Mpp, int dv, int gs) {
    const Model& m = Mdl();
    if (dv == kDirOff) return 1.;
    const int v = dv >= kV3 ? 2 : dv - 1;
    if (Mpp <= 2 * kMpi || Mpp >= kMetaGen - 2 * kMe) return 0.;
    const double tp = Interp(m.Tpp[v][gs], m.m0, m.m1, Mpp) / m.ZT[v][gs];
    const double gp = m.h[BinOf(m.mE, Mpp)] * Interp(m.Glam, m.m0, m.m1, Mpp);
    return gp > 0 ? m.renPP[v][gs] * tp / gp : 0.;
}
inline double WBump(double Mpp, int bv) {
    if (bv == kBumpOff) return 1.;
    return DalFbar(Mpp, bv) / Mdl().normB[bv];
}
inline double WConv(double Mpp, int gs) {
    const double ps = PStar(Mpp, kMpi), k = K3(kMetaHG, Mpp, 0.);
    return ps * ps * k * k * FV2t(Mpp, gs) / Mdl().normC[gs];
}
// HGeant probability that a pion pair of one eta with this M(pipi) comes from pi+pi-pi0
inline double P3(double Mpp) {
    const Model& m = Mdl();
    const double a = kHGppp0 * Rho3(Mpp) / m.rho3N, b = kHGppg * RhoG(Mpp) / m.rhoGN;
    return (a + b) > 0 ? a / (a + b) : 0.;
}
inline double WPcls2(double Mpp, int bv, int conv, int gs) {
    const double p3 = P3(Mpp);
    const double wb = (bv && p3 > 0) ? WBump(Mpp, bv) : 1.;
    const double wc = conv ? WConv(Mpp, gs) : 1.;
    return p3 * wb + (1 - p3) * wc;
}

// ---------------------------------------------------------------- true angles of class-5 rows (variant V3)
struct AngKey {
    int grp; float a, b, c, d;
    bool operator==(const AngKey& o) const { return grp == o.grp && a == o.a && b == o.b && c == o.c && d == o.d; }
};
struct AngHash {
    size_t operator()(const AngKey& k) const {
        unsigned int u[4];
        std::memcpy(&u[0], &k.a, 4); std::memcpy(&u[1], &k.b, 4); std::memcpy(&u[2], &k.c, 4); std::memcpy(&u[3], &k.d, 4);
        unsigned long long h = 1469598103934665603ULL ^ (unsigned long long)(unsigned)k.grp;
        for (int i = 0; i < 4; ++i) { h ^= u[i]; h *= 1099511628211ULL; }
        return (size_t)h;
    }
};
struct AngVal { double ce, cp, c2p; };
inline std::unordered_map<AngKey, AngVal, AngHash>& AngMap() { static std::unordered_map<AngKey, AngVal, AngHash> m; return m; }
inline long long& AngMiss() { static long long n = 0; return n; }
inline long long& AngHit() { static long long n = 0; return n; }
inline bool& AngHavePhi() { static bool b = false; return b; }
inline bool& AngTried() { static bool b = false; return b; }
inline bool LoadAngles(const char* file = kRwAngFile) {
    AngTried() = true;
    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie()) { printf("simphys::LoadAngles: cannot open %s\n", file); return false; }
    TTree* t = (TTree*)f->Get("ang");
    if (!t) { printf("simphys::LoadAngles: tree 'ang' missing\n"); return false; }
    Int_t grp; Float_t a, b, c, d; Double_t ce, cp, c2p = 0;
    t->SetBranchAddress("grp", &grp); t->SetBranchAddress("ep_sim_p", &a); t->SetBranchAddress("em_sim_p", &b);
    t->SetBranchAddress("pip_sim_p", &c); t->SetBranchAddress("pim_sim_p", &d);
    t->SetBranchAddress("ce", &ce); t->SetBranchAddress("cp", &cp);
    AngHavePhi() = t->GetBranch("c2phi") != nullptr;
    if (AngHavePhi()) t->SetBranchAddress("c2phi", &c2p);
    AngMap().clear();
    AngMap().reserve((size_t)t->GetEntries() * 2);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); AngMap()[AngKey{grp, a, b, c, d}] = AngVal{ce, cp, c2p}; }
    printf("simphys::LoadAngles: %zu class-5 rows with true angles%s\n", AngMap().size(), AngHavePhi() ? " and cos(2 phi)" : " (no cos(2 phi): V4 = V3)");
    f->Close(); delete f;
    return true;
}
inline bool FindAngles(const Row& r, double& ce, double& cp, double& c2p) {
    if (AngMap().empty() && !AngTried()) {   // V3 is the main variant: load the angles on first use
        printf("simphys: V3/V4 weight requested before LoadAngles(): loading %s\n", kRwAngFile);
        LoadAngles();
    }
    if (AngMap().empty()) {
        static bool warned = false;
        if (!warned) { printf("simphys: WARNING no true angles available: V3/V4 fall back to the V2 shape\n"); warned = true; }
        ++AngMiss();
        return false;
    }
    auto it = AngMap().find(AngKey{r.grp, r.ep_sim_p, r.em_sim_p, r.pip_sim_p, r.pim_sim_p});
    if (it == AngMap().end()) { ++AngMiss(); return false; }
    ++AngHit(); ce = it->second.ce; cp = it->second.cp; c2p = it->second.c2p;
    return true;
}

// ---------------------------------------------------------------- the factor of one row
inline double RowFactor(const Row& r, int code) {
    if (!code) return 1.;
    if (r.cls == 6) return 1.;                       // X17: unweighted
    const int dv = CodeDir(code), bv = CodeBump(code), cv = CodeConv(code), bk = CodeBkg(code), gs = CodeGS(code);
    double f = 1.;
    if (r.set == 3) {
        if (r.cls == 5) {
            if (dv) {
                double ce = 0, cp = 0, c2p = 0;
                const bool ha = (dv >= kV3) && FindAngles(r, ce, cp, c2p);
                f *= W5(r.mee_true, r.mpipi_true, dv, gs, ce, cp, ha, c2p);
            }
            return f;
        }
        if (dv && bk) {
            if (r.lcls == 10 && r.mee_true > 0) f *= Wee(r.mee_true, dv, gs);
            if (r.pcls == 9 && r.mpipi_true > 0) f *= Wpp(r.mpipi_true, dv, gs);
        }
    }
    if (r.cls == 2 || r.cls == 3) { if (bv) f *= WBump(r.mpipi_true, bv); }
    else if (r.cls == 4) { if (cv) f *= WConv(r.mpipi_true, gs); }
    else if (r.pcls == 2 && bk && (bv || cv) && r.mpipi_true > 0) f *= WPcls2(r.mpipi_true, bv, cv, gs);
    return f;
}
inline double RowFactorHook(const Row& r, int code) { return RowFactor(r, code); }

// ---------------------------------------------------------------- initialisation
// Gauss-Legendre nodes on [-1, 1]
inline void GL(int n, std::vector<double>& x, std::vector<double>& w) {
    x.assign(n, 0.); w.assign(n, 0.);
    for (int i = 0; i < n; ++i) {
        double z = std::cos(M_PI * (i + 0.75) / (n + 0.5)), z1, pp;
        do {
            double p1 = 1, p2 = 0;
            for (int j = 0; j < n; ++j) { const double p3 = p2; p2 = p1; p1 = ((2 * j + 1) * z * p2 - j * p3) / (j + 1); }
            pp = n * (z * p1 - p2) / (z * z - 1);
            z1 = z; z = z1 - p1 / pp;
        } while (std::fabs(z - z1) > 1e-14);
        x[i] = z; w[i] = 2 / ((1 - z * z) * pp * pp);
    }
}
// read G and the stored constants; build the fine tables and the T normalisations
inline bool Init(const char* file = kRwGenFile, bool verbose = true) {
    Model& m = Mdl();
    m.ok = false;
    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie()) { printf("simphys::Init: cannot open %s\n", file); return false; }
    TH1D* hg = (TH1D*)f->Get("g_q");
    TH1D* hh = (TH1D*)f->Get("h_m");
    TH1D* hc = (TH1D*)f->Get("consts");
    if (!hg || !hh || !hc) { printf("simphys::Init: g_q, h_m or consts missing in %s\n", file); return false; }
    m.qE.clear(); m.g.clear(); m.mE.clear(); m.h.clear();
    for (int i = 1; i <= hg->GetNbinsX() + 1; ++i) m.qE.push_back(hg->GetXaxis()->GetBinLowEdge(i));
    for (int i = 1; i <= hg->GetNbinsX(); ++i) m.g.push_back(hg->GetBinContent(i));
    for (int i = 1; i <= hh->GetNbinsX() + 1; ++i) m.mE.push_back(hh->GetXaxis()->GetBinLowEdge(i));
    for (int i = 1; i <= hh->GetNbinsX(); ++i) m.h.push_back(hh->GetBinContent(i));
    // consts: 1-8 ren5[v][fv] (v 0..3), 9-14 renEE, 15-20 renPP, 21 bE, 22 bPi, 23 normB KLOE, 24 normB slides
    auto C = [&](int i) { return hc->GetBinContent(i); };
    for (int v = 0; v < 4; ++v) for (int fv = 0; fv < 2; ++fv) m.ren5[v][fv] = C(1 + 2 * v + fv) > 0 ? C(1 + 2 * v + fv) : 1.;
    for (int v = 0; v < 3; ++v) for (int fv = 0; fv < 2; ++fv) {
        m.renEE[v][fv] = C(9 + 2 * v + fv) > 0 ? C(9 + 2 * v + fv) : 1.;
        m.renPP[v][fv] = C(15 + 2 * v + fv) > 0 ? C(15 + 2 * v + fv) : 1.;
    }
    m.bE = C(21); m.bPi = C(22);
    m.normB[kBumpKLOE] = C(23) > 0 ? C(23) : 1.; m.normB[kBumpSlides] = C(24) > 0 ? C(24) : 1.;
    f->Close(); delete f;

    // |F_V|^2 table
    m.fvN = 6001; m.fvM0 = 2 * kMpi; m.fvM1 = 0.56;
    for (int gs = 0; gs < 2; ++gs) {
        m.fvTab[gs].assign(m.fvN, 0.);
        for (int i = 0; i < m.fvN; ++i) { const double x = m.fvM0 + (m.fvM1 - m.fvM0) * i / (m.fvN - 1); m.fvTab[gs][i] = FV2(x * x, gs); }
    }
    // threshold bin average of KW
    {
        const double a = m.qE[0], b = m.qE[1];
        const int n = 20000; double s = 0;
        const double smax = std::sqrt(b - a);
        for (int i = 0; i < n; ++i) { const double sv = smax * (i + 0.5) / n; s += KW(a + sv * sv) * 2 * sv * (smax / n); }
        m.KWbar0 = s / (b - a);
    }
    std::vector<double> gx, gw, hx, hw;
    GL(64, gx, gw);
    GL(16, hx, hw);
    // fine q grid (log): Hlam, Iee
    m.nQ = 3001; m.lq0 = std::log(2 * kMe); m.lq1 = std::log(QMax());
    m.Hlam.assign(m.nQ, 0.);
    for (int v = 0; v < 3; ++v) for (int gs = 0; gs < 2; ++gs) m.Iee[v][gs].assign(m.nQ, 0.);
    for (int i = 0; i < m.nQ; ++i) {
        const double q = std::exp(m.lq0 + (m.lq1 - m.lq0) * i / (m.nQ - 1));
        const double lo = 2 * kMpi, hi = kMetaGen - q;
        if (hi <= lo) continue;
        // H(q): integrate the piecewise h exactly bin by bin (GL inside each bin part)
        double H = 0;
        for (size_t b = 0; b + 1 < m.mE.size(); ++b) {
            const double a0 = std::max(lo, m.mE[b]), a1 = std::min(hi, m.mE[b + 1]);
            if (a1 <= a0) continue;
            double s = 0;
            for (int k = 0; k < 16; ++k) {
                const double xx = 0.5 * (a0 + a1) + 0.5 * (a1 - a0) * hx[k];
                const double l = Lam(kMetaGen * kMetaGen, xx * xx, q * q);
                s += 0.5 * (a1 - a0) * hw[k] * (l > 0 ? std::sqrt(l) : 0.);
            }
            H += m.h[b] * s;
        }
        m.Hlam[i] = H;
        double s3[3][2] = {{0, 0}, {0, 0}, {0, 0}};
        for (int k = 0; k < 64; ++k) {
            const double xx = 0.5 * (lo + hi) + 0.5 * (hi - lo) * gx[k];
            const double ww = 0.5 * (hi - lo) * gw[k];
            const double r0 = Rmass(q, xx, 0, 0);
            for (int gs = 0; gs < 2; ++gs) {
                const double f2 = FV2t(xx, gs);
                s3[0][gs] += ww * r0;
                s3[1][gs] += ww * r0 * f2;
                s3[2][gs] += ww * r0 * f2 * FVMD2(q * q);
            }
        }
        for (int v = 0; v < 3; ++v) for (int gs = 0; gs < 2; ++gs) m.Iee[v][gs][i] = s3[v][gs];
    }
    // normalisations ZT: int dq KWeff(q) int dM R (q = 2 m_e + s^2)
    {
        const int n = 40000;
        const double a = 2 * kMe, smax = std::sqrt(QMax() - a);
        double Z[3][2] = {{0, 0}, {0, 0}, {0, 0}};
        for (int i = 0; i < n; ++i) {
            const double sv = smax * (i + 0.5) / n, q = a + sv * sv, dq = 2 * sv * smax / n;
            const double lq = std::log(q);
            for (int v = 0; v < 3; ++v) for (int gs = 0; gs < 2; ++gs) {
                Z[v][gs] += KWeff(q) * Interp(m.Iee[v][gs], m.lq0, m.lq1, lq) * dq;
            }
        }
        for (int v = 0; v < 3; ++v) for (int gs = 0; gs < 2; ++gs) m.ZT[v][gs] = Z[v][gs];
    }
    // fine M grid: Glam = int g lambda^{1/2} dq, Tpp = int T dq (unnormalised)
    m.nM = 1201; m.m0 = 2 * kMpi; m.m1 = kMetaGen - 2 * kMe;
    m.Glam.assign(m.nM, 0.);
    for (int v = 0; v < 3; ++v) for (int gs = 0; gs < 2; ++gs) m.Tpp[v][gs].assign(m.nM, 0.);
    for (int j = 0; j < m.nM; ++j) {
        const double Mpp = m.m0 + (m.m1 - m.m0) * j / (m.nM - 1);
        const double qhi = kMetaGen - Mpp, a = 2 * kMe;
        if (qhi <= a) continue;
        const int n = 3000;
        const double smax = std::sqrt(qhi - a);
        double G = 0, T[3][2] = {{0, 0}, {0, 0}, {0, 0}};
        const double f2[2] = {FV2t(Mpp, 0), FV2t(Mpp, 1)};
        for (int i = 0; i < n; ++i) {
            const double sv = smax * (i + 0.5) / n, q = a + sv * sv, dq = 2 * sv * smax / n;
            const double l = Lam(kMetaGen * kMetaGen, Mpp * Mpp, q * q);
            if (l <= 0) continue;
            G += Gq(q) * std::sqrt(l) * dq;
            const double r0 = KWeff(q) * Rmass(q, Mpp, 0, 0) * dq;
            for (int gs = 0; gs < 2; ++gs) {
                T[0][gs] += r0; T[1][gs] += r0 * f2[gs]; T[2][gs] += r0 * f2[gs] * FVMD2(q * q);
            }
        }
        m.Glam[j] = G;
        for (int v = 0; v < 3; ++v) for (int gs = 0; gs < 2; ++gs) m.Tpp[v][gs][j] = T[v][gs];
    }
    // conversion normalisation <|F_V|^2 p*^2 k^2>_PS and the phase-space normalisations of Rho3, RhoG
    {
        const int n = 20000;
        double sPS = 0, s3 = 0, sC[2] = {0, 0};
        for (int i = 0; i < n; ++i) {
            const double x = 2 * kMpi + (kMetaHG - 2 * kMpi) * (i + 0.5) / n, dx = (kMetaHG - 2 * kMpi) / n;
            const double rg = RhoG(x);
            sPS += rg * dx; s3 += Rho3(x) * dx;
            const double ps = PStar(x, kMpi), k = K3(kMetaHG, x, 0.);
            for (int gs = 0; gs < 2; ++gs) sC[gs] += rg * ps * ps * k * k * FV2t(x, gs) * dx;
        }
        m.rhoGN = sPS; m.rho3N = s3;
        for (int gs = 0; gs < 2; ++gs) m.normC[gs] = sC[gs] / sPS;
    }
    m.ok = true;
    if (verbose)
        printf("simphys::Init: G with %zu q bins x %zu M bins; ZT(V0) %.4e, ZT(V2, GP) %.4e; <f>_PS KLOE %.5f, slides %.5f; "
               "ren5(V2) %.5f; bE %.4f, bPi %.4f\n",
               m.g.size(), m.h.size(), m.ZT[0][0], m.ZT[2][0], m.normB[1], m.normB[2], m.ren5[2][0], m.bE, m.bPi);
    return true;
}
// connect to SimEngine.h (SimOpt::phys)
inline void Install() { simeng::PhysHook() = &RowFactorHook; }

}  // namespace simphys

#endif  // SIM_ETA_PHYSICS_H
