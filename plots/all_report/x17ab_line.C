// ============================================================================
// x17ab_line.C  --  X17 line shapes and efficiencies for every hypothetical mass (SPEC report_all, WP9 X-a)
// ============================================================================
// HADES p+p, T = 4.5 GeV, eta -> pi+ pi- X17, X17 -> e+ e-.  SIMULATION ONLY (eta set, owned rows, mode A).
// Selection S0 (pre + both pions in the vertex fit), simulated M4' = s [m4_true + r (m4 - m4_true)] (s = 0.99615,
// r = 1) inside the X2D M4 range 0.5225-0.6100 GeV/c2 (data scale).
//
// What is done, in plain words:
//  1. The only simulated X17 sample (class 6) has m_X = 17.000 MeV/c2. Its reconstructed M(ee) line (0.25 MeV bins),
//     its M4' shape and its efficiency eps_X(17) = sigma_vis per unit BR / sigma(pp -> eta X) are stored.
//  2. For any other mass m the line is built from the direct decays eta -> pi+pi-e+e- (class 5) whose TRUE M(ee) lies
//     within +-0.5 MeV/c2 of m. Each such row is
//       (a) shifted in reconstructed M(ee) by (m - q_true), so that the slice width does not broaden the line;
//       (b) reweighted from the generator kinematics to the X17 kinematics at fixed q:
//           M(pipi): three-body phase space of eta -> pi+ pi- X (density p*(M) lambda^1/2(m_eta^2, M^2, q^2)) instead
//                    of the generator density h(M) lambda^1/2 (sim_eta_physics.h, G = g(q) h(M) lambda^1/2);
//           e+ decay angle in the e+e- frame: isotropic (the generator is flat already, factor 1/a_e(c));
//           pi+ decay angle in the pi+pi- frame: isotropic instead of the generator's sin^2 (weight 2/(3(1-c^2)),
//                    capped at |c| = 0.97, the capped mean is divided out).
//     The weights average to 1 over the generated decays at fixed q, so the number of generated decays of the slice is
//     the generator density integrated over the slice: N_gen(slice) = sigma_eta BR_dir L_sim x int_slice g(q) H(q) dq
//     (H(q) = int h lambda^1/2 dM). The slice efficiency eps_slice(m) = sum of the accepted weights / N_gen(slice).
//  3. Variants of the reweighting (nuisance "eps_X/eps_5" of the spec):
//       DL   M(pipi) "direct-like" (the V2 direct-decay density at fixed q, w5(V2)/Wee(V2)), decay angles isotropic;
//       C2   phase-space M(pipi), e+ decay angle 1 + cos^2 (transverse virtual photon: V3 angular factor);
//       NP   no pion-angle reweighting (generator sin^2 kept).
//  4. Validation on C6 at 17 MeV/c2: line mean, RMS, double-Gaussian fit, chi2 of the two shapes, eps_slice(17)/eps_X(17).
//  5. The model uses eps_X(m) = eps_X(17, C6) x eps_slice(m) / eps_slice(17) (anchored on C6) and the line shape of the
//     slice (C6 itself at 17.0); R_X(m) = eps_X(m) / eps_5, eps_5 = C5 (V3) / (sigma_eta BR_dir), same selection and range.
//
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/x17ab/x17_lines.root (histograms per mass and variant, TTree "lines"),
//         docs/latex/tables_all/x17ab_lines.txt, docs/latex/figs_all/x17ab_line_valid.pdf/.png, x17ab_line_eff.pdf/.png
// Run (from /home/przygoda/HADES/FAT/FAT; about 3 GB, 2 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17ab", kTRUE)' 'plots/all_report/x17ab_line.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "defs.h"
#include "sim_eta_physics.h"

namespace x17line {
using namespace alldefs;

const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/x17ab";
const char* const kFigDir = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
const char* const kTabDir = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all";
const double kSigEtaNb = kSigEtaSMASH * 1e6;   // nb
const double kM4lo = kX2Dxlo, kM4hi = kX2Dxhi;  // data-scale M4' range of the model
// mass grid [MeV/c2]
const double kMlo = 6.0, kMhi = 50.0, kMstep = 0.5;
const int kNM = (int)std::lround((kMhi - kMlo) / kMstep) + 1;
const double kHalf = 0.5;                       // slice half width [MeV/c2]
// variants
enum { kVnomX = 0, kVDL, kVC2, kVNP, kNVX };
const char* kVname[kNVX] = {"nom", "DL", "C2", "NP"};
const char* kVtext[kNVX] = {"phase-space M(pipi), isotropic e+ and pi+ angles",
                            "direct-like M(pipi) (V2 at fixed q), isotropic angles",
                            "phase-space M(pipi), e+ angle 1+cos^2 (V3 factor)",
                            "phase-space M(pipi), isotropic e+, generator pi+ angle (sin^2)"};
// fine M(ee) histogram
const int kNee = 400;
const double kEeLo = 0., kEeHi = 100.;   // MeV/c2, 0.25 MeV bins
const double kPiCap = 0.97;

// phase-space conditional density in M at fixed q (unnormalised) and its normalisation
inline double PSm(double M, double q) {
    const double me = simphys::kMetaGen;
    const double l = simphys::Lam(me * me, M * M, q * q);
    if (l <= 0 || M <= 2 * simphys::kMpi) return 0.;
    return simphys::PStar(M, simphys::kMpi) * std::sqrt(l);
}
inline double ZPS(double q) {   // int PSm dM
    const double lo = 2 * simphys::kMpi, hi = simphys::kMetaGen - q;
    if (hi <= lo) return 0.;
    const int n = 2000;
    double s = 0;
    for (int i = 0; i < n; ++i) { const double M = lo + (hi - lo) * (i + 0.5) / n; s += PSm(M, q); }
    return s * (hi - lo) / n;
}
inline double Hq(double q) { const simphys::Model& m = simphys::Mdl(); return simphys::Interp(m.Hlam, m.lq0, m.lq1, std::log(q)); }
// generator q-marginal density g(q) H(q) integrated over [a, b] (GeV)
inline double PgenQ(double a, double b) {
    const int n = 400;
    double s = 0;
    for (int i = 0; i < n; ++i) { const double q = a + (b - a) * (i + 0.5) / n; s += simphys::Gq(q) * Hq(q); }
    return s * (b - a) / n;
}
// mean of the capped pion weight under sin^2: int_{|c|<cap} 0.5 dc + int_{|c|>cap} 0.75(1-c^2) wcap dc
inline double PiCapMean() {
    const double wcap = 2. / (3. * (1 - kPiCap * kPiCap));
    const double inner = kPiCap;   // 2 x 0.5 x cap
    // 2 x int_cap^1 0.75 (1 - c^2) dc = 1.5 [ (1 - cap) - (1 - cap^3)/3 ]
    const double outer = 1.5 * ((1 - kPiCap) - (1 - kPiCap * kPiCap * kPiCap) / 3.) * wcap;
    return inner + outer;
}
inline double WPi(double cp) {
    const double c = std::min(std::fabs(cp), kPiCap);
    return 2. / (3. * (1 - c * c));
}

struct Acc {
    TH1D* hee = nullptr;   // reco M(ee), shifted (MeV)
    TH1D* hm4 = nullptr;   // M4' (GeV)
    double sw = 0, sw2 = 0;
    long long n = 0;
    void Init(const char* name) {
        hee = new TH1D(Form("ee_%s", name), "", kNee, kEeLo, kEeHi); hee->Sumw2(); hee->SetDirectory(nullptr);
        hm4 = new TH1D(Form("m4_%s", name), "", kX2Dnx, kM4lo, kM4hi); hm4->Sumw2(); hm4->SetDirectory(nullptr);
    }
    void Fill(double mee, double m4, double w) { hee->Fill(mee, w); hm4->Fill(m4, w); sw += w; sw2 += w * w; ++n; }
};

// double-Gaussian fit of a line (unit area); returns parameters mean, s1, s2, f1 and chi2/ndf
struct DG { double mu = 0, s1 = 0, s2 = 0, f1 = 1, chi2 = 0; int ndf = 0; double mean = 0, rms = 0; };
DG FitDG(TH1D* h, double m) {
    DG r;
    if (!h || h->Integral() <= 0) return r;
    TH1D* c = (TH1D*)h->Clone(Form("%s_fit", h->GetName()));
    c->SetDirectory(nullptr);
    c->Scale(1. / c->Integral() / c->GetBinWidth(1));
    r.mean = h->GetMean(); r.rms = h->GetStdDev();
    TF1 f("dg", "[0]*([3]*TMath::Gaus(x,[1],[2],1)+(1-[3])*TMath::Gaus(x,[1],[4],1))", std::max(0.5, m - 15.), m + 15.);
    f.SetParameters(1., m, 0.6 * r.rms, 0.7, 1.6 * r.rms);
    f.SetParLimits(3, 0.05, 1.0);
    f.SetParLimits(2, 0.3, 30.); f.SetParLimits(4, 0.3, 40.);
    f.FixParameter(0, 1.);
    c->Fit(&f, "QNR0");
    r.mu = f.GetParameter(1); r.s1 = f.GetParameter(2); r.f1 = f.GetParameter(3); r.s2 = f.GetParameter(4);
    if (r.s1 > r.s2) { std::swap(r.s1, r.s2); r.f1 = 1 - r.f1; }
    r.chi2 = f.GetChisquare(); r.ndf = f.GetNDF();
    delete c;
    return r;
}
// chi2 between two unit-area shapes (fine bins merged by 4 = 1 MeV), using their errors, in [lo, hi] MeV
double ShapeChi2(TH1D* a, TH1D* b, double lo, double hi, int& ndf) {
    TH1D* x = (TH1D*)a->Clone("sca"); TH1D* y = (TH1D*)b->Clone("scb");
    x->SetDirectory(nullptr); y->SetDirectory(nullptr);
    x->Rebin(4); y->Rebin(4);
    x->Scale(1. / x->Integral()); y->Scale(1. / y->Integral());
    double c2 = 0; ndf = 0;
    for (int i = 1; i <= x->GetNbinsX(); ++i) {
        const double ce = x->GetBinCenter(i);
        if (ce < lo || ce > hi) continue;
        if (x->GetBinContent(i) <= 0 || y->GetBinContent(i) <= 0) continue;   // empty bins of the weighted slice have no error
        const double e2 = x->GetBinError(i) * x->GetBinError(i) + y->GetBinError(i) * y->GetBinError(i);
        if (e2 <= 0) continue;
        const double d = x->GetBinContent(i) - y->GetBinContent(i);
        c2 += d * d / e2; ++ndf;
    }
    delete x; delete y;
    return c2;
}

// mean and RMS inside [m - 12, m + 12] MeV/c2 (robust against the far tails)
void TruncMoments(TH1D* h, double m, double& mean, double& rms) {
    double s0 = 0, s1 = 0, s2 = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        const double x = h->GetBinCenter(i), w = h->GetBinContent(i);
        if (x < m - 12 || x > m + 12) continue;
        s0 += w; s1 += w * x; s2 += w * x * x;
    }
    mean = s0 > 0 ? s1 / s0 : 0; rms = s0 > 0 ? std::sqrt(std::max(0., s2 / s0 - mean * mean)) : 0;
}

}  // namespace x17line

void x17ab_line() {
    using namespace x17line;
    gSystem->mkdir(kOutDir, kTRUE);
    gStyle->SetOptStat(0);
    simphys::Init();
    simphys::LoadAngles();
    const simphys::Model& md = simphys::Mdl();
    printf("x17ab_line: generator e+ angle b = %.4f, pion b = %.1f (-1 = sin^2); capped pion-weight mean %.5f\n", md.bE, md.bPi, PiCapMean());
    // normalisation check of the generator q-marginal
    const double ptot = PgenQ(2 * simphys::kMe * 1.000001, simphys::QMax());
    printf("x17ab_line: int g(q) H(q) dq over the full range = %.5f (expected 1)\n", ptot);

    simeng::SimData* D = new simeng::SimData();
    if (!simeng::Load(*D, simeng::kCompact, 1 << simeng::kEtaUS, 0.50, 0.65, 1)) return;
    const std::vector<simeng::Row>& v = D->t[simeng::kEtaUS];
    printf("x17ab_line: %zu eta_us rows in raw M4 0.50-0.65\n", v.size());

    // accumulators: C6 (X17, 17 MeV), C5 V3 (D template), per mass x variant slices
    Acc c6, c5;
    c6.Init("C6"); c5.Init("C5V3");
    std::vector<std::vector<Acc>> sl(kNM, std::vector<Acc>(kNVX));
    for (int k = 0; k < kNM; ++k) for (int x = 0; x < kNVX; ++x) sl[k][x].Init(Form("m%03d_%s", (int)std::lround(10 * (kMlo + k * kMstep)), kVname[x]));
    std::map<long long, double> zcache;   // ZPS(q) cache on a 1 keV grid
    auto ZPSc = [&](double q) {
        const long long key = std::llround(q * 1e6);
        auto it = zcache.find(key);
        if (it != zcache.end()) return it->second;
        const double z = ZPS(key * 1e-6);
        zcache[key] = z;
        return z;
    };
    const double piMean = PiCapMean();
    long long nC5 = 0, nC6 = 0, nNoAng = 0;
    for (const auto& r : v) {
        if (r.ownA != 1 || r.tmpl != 1) continue;
        if (r.cls != 5 && r.cls != 6) continue;
        if (!PassSel(kSelS0, r, r.cor_ok == 1, DzPi(r.pip_z, r.pim_z, r.vz), r.mm4)) continue;
        const double m4p = M4sim(r.m4, r.m4_true);
        if (!(m4p >= kM4lo && m4p < kM4hi)) continue;
        const double meeR = 1000. * r.mee;   // MeV
        if (r.cls == 6) { c6.Fill(meeR, m4p, 2.); ++nC6; continue; }
        ++nC5;
        const double w3 = r.w * simphys::RowFactor(r, kPhysNom);
        c5.Fill(meeR, m4p, w3);
        const double q = r.mee_true, M = r.mpipi_true;   // GeV
        if (!(q > 0) || !(M > 0)) continue;
        const double qM = 1000. * q;
        const int k0 = std::max(0, (int)std::ceil((qM - kHalf - kMlo) / kMstep - 1e-9));
        const int k1 = std::min(kNM - 1, (int)std::floor((qM + kHalf - kMlo) / kMstep + 1e-9));
        if (k1 < k0) continue;
        // weights at fixed q
        double ce = 0, cp = 0, c2p = 0;
        const bool ha = simphys::FindAngles(r, ce, cp, c2p);
        if (!ha) { ++nNoAng; }
        const double hm = md.h[simphys::BinOf(md.mE, M)];
        const double Hv = Hq(q), Zv = ZPSc(q);
        const double wps = (hm > 0 && Zv > 0) ? simphys::PStar(M, simphys::kMpi) * Hv / (hm * Zv) : 0.;
        const double ae = (1 + md.bE * ce * ce) / (2 + 2 * md.bE / 3);
        const double wE = ha ? 0.5 / ae : 1.;
        const double wEc2 = ha ? simphys::AngE(ce, q) / ae : 1.;
        const double wP = ha ? WPi(cp) / piMean : 1.;
        // direct-like M(pipi) at fixed q: V2 weight divided by its q-marginal (V2 = no angle factors)
        const double w5v2 = simphys::W5(q, M, simphys::kV2, 0);
        const double wee2 = simphys::Wee(q, simphys::kV2, 0);
        const double wdl = wee2 > 0 ? w5v2 / wee2 : 0.;
        const double wx[kNVX] = {wps * wE * wP, wdl * wE * wP, wps * wEc2 * wP, wps * wE};
        for (int k = k0; k <= k1; ++k) {
            const double m = kMlo + k * kMstep;
            if (std::fabs(qM - m) >= kHalf) continue;
            const double meeS = meeR + (m - qM);
            for (int x = 0; x < kNVX; ++x) sl[k][x].Fill(meeS, m4p, r.w * wx[x]);
        }
    }
    printf("x17ab_line: accepted C5 rows %lld (no angles %lld), C6 rows %lld\n", nC5, nNoAng, nC6);
    delete D;

    // efficiencies
    const double Lsim = kLsim;
    const double sigC6 = c6.sw / Lsim;                      // nb per unit BR (S0, X2D range, all M(ee))
    const double epsX17 = sigC6 / kSigEtaNb;
    const double sigC5 = c5.sw / Lsim;                      // nb (V3)
    const double eps5 = sigC5 / (kSigEtaNb * kBRdir);
    printf("x17ab_line: C6 %.3f nb per unit BR -> eps_X(17) = %.5e; C5(V3) %.4f nb -> eps_5 = %.5e; R_X(17, C6) = %.4f\n", sigC6, epsX17, sigC5, eps5, epsX17 / eps5);
    // slices
    std::vector<double> epsS[kNVX], epsSe[kNVX], pgen(kNM);
    for (int k = 0; k < kNM; ++k) {
        const double m = kMlo + k * kMstep;
        pgen[k] = PgenQ(1e-3 * (m - kHalf), 1e-3 * (m + kHalf));
        for (int x = 0; x < kNVX; ++x) {
            const double sgen = kSigEtaNb * kBRdir * pgen[k] * Lsim;   // generated rows x weight per unit generated decay
            // r.w of a class-5 row = (2 BR_dir per forced decay); accepted sum of r.w x wx / L_sim = sigma_vis [nb]
            const double e = sl[k][x].sw / Lsim / (kSigEtaNb * kBRdir * pgen[k]);
            const double ee = std::sqrt(sl[k][x].sw2) / Lsim / (kSigEtaNb * kBRdir * pgen[k]);
            (void)sgen;
            epsS[x].push_back(e); epsSe[x].push_back(ee);
        }
    }
    const int k17 = (int)std::lround((17.0 - kMlo) / kMstep);
    printf("x17ab_line: eps_slice(17) nom = %.5e +- %.2e; ratio to C6 = %.4f +- %.4f\n", epsS[0][k17], epsSe[0][k17], epsS[0][k17] / epsX17, epsSe[0][k17] / epsX17);

    // fits and validation
    DG dg6 = FitDG(c6.hee, 17.);
    std::vector<DG> dgs(kNM);
    for (int k = 0; k < kNM; ++k) dgs[k] = FitDG(sl[k][0].hee, kMlo + k * kMstep);
    int ndfS = 0;
    const double chiS = ShapeChi2(c6.hee, sl[k17][0].hee, 5., 32., ndfS);
    int ndfM = 0;
    // M4 shape comparison C6 vs C5 (V3) vs slice
    auto M4chi = [&](TH1D* a, TH1D* b, int& nd) {
        double c2 = 0; nd = 0;
        const double ia = a->Integral(), ib = b->Integral();
        for (int i = 1; i <= a->GetNbinsX(); ++i) {
            if (a->GetBinContent(i) <= 0 || b->GetBinContent(i) <= 0) continue;
            const double ea = a->GetBinError(i) / ia, eb = b->GetBinError(i) / ib;
            const double e2 = ea * ea + eb * eb;
            if (e2 <= 0) continue;
            const double d = a->GetBinContent(i) / ia - b->GetBinContent(i) / ib;
            c2 += d * d / e2; ++nd;
        }
        return c2;
    };
    const double chiM4a = M4chi(c6.hm4, c5.hm4, ndfM);
    int ndfM2 = 0;
    const double chiM4b = M4chi(c6.hm4, sl[k17][0].hm4, ndfM2);
    printf("x17ab_line: C6 line: mean %.3f rms %.3f; DG mu %.3f s1 %.3f s2 %.3f f1 %.3f\n", dg6.mean, dg6.rms, dg6.mu, dg6.s1, dg6.s2, dg6.f1);
    printf("x17ab_line: slice(17) line: mean %.3f rms %.3f; DG mu %.3f s1 %.3f s2 %.3f f1 %.3f; shape chi2 (5-32 MeV, 1 MeV bins, non-empty) %.1f / %d\n",
           dgs[k17].mean, dgs[k17].rms, dgs[k17].mu, dgs[k17].s1, dgs[k17].s2, dgs[k17].f1, chiS, ndfS);
    printf("x17ab_line: M4' shape chi2: C6 vs C5(V3) %.1f / %d; C6 vs slice(17) %.1f / %d\n", chiM4a, ndfM, chiM4b, ndfM2);

    // smooth mass dependence: truncated mean and RMS of the slice lines and the slice efficiencies, pol2 in m over 8-45 MeV
    double tm6 = 0, tr6 = 0;
    TruncMoments(c6.hee, 17., tm6, tr6);
    TGraphErrors gMean, gRms, gEps[kNVX];
    std::vector<double> tmS(kNM), trS(kNM);
    for (int k = 0; k < kNM; ++k) {
        const double m = kMlo + k * kMstep;
        TruncMoments(sl[k][0].hee, m, tmS[k], trS[k]);
        if (m < 8 - 1e-9 || m > 45 + 1e-9) continue;
        const double ne = std::max(10., sl[k][0].hee->GetEffectiveEntries());
        gMean.SetPoint(gMean.GetN(), m, tmS[k] - m); gMean.SetPointError(gMean.GetN() - 1, 0, trS[k] / std::sqrt(ne));
        gRms.SetPoint(gRms.GetN(), m, trS[k]); gRms.SetPointError(gRms.GetN() - 1, 0, trS[k] / std::sqrt(2 * ne));
        for (int x = 0; x < kNVX; ++x) { gEps[x].SetPoint(gEps[x].GetN(), m, epsS[x][k]); gEps[x].SetPointError(gEps[x].GetN() - 1, 0, epsSe[x][k]); }
    }
    TF1 fMean("fMean", "pol2", 8, 45), fRms("fRms", "pol2", 8, 45);
    gMean.Fit(&fMean, "QN0"); gRms.Fit(&fRms, "QN0");
    TF1* fEps[kNVX];
    for (int x = 0; x < kNVX; ++x) { fEps[x] = new TF1(Form("fEps_%s", kVname[x]), "pol2", 8, 45); gEps[x].Fit(fEps[x], "QN0"); }
    printf("x17ab_line: smooth: mean-m = %.3f%+.4f m%+.6f m^2 (chi2 %.1f/%d), rms = %.3f%+.4f m%+.6f m^2 (chi2 %.1f/%d); C6 trunc mean %.3f rms %.3f\n",
           fMean.GetParameter(0), fMean.GetParameter(1), fMean.GetParameter(2), fMean.GetChisquare(), fMean.GetNDF(), fRms.GetParameter(0),
           fRms.GetParameter(1), fRms.GetParameter(2), fRms.GetChisquare(), fRms.GetNDF(), tm6, tr6);
    for (int x = 0; x < kNVX; ++x)
        printf("x17ab_line: smooth eps %s: chi2 %.1f/%d; R_X(10, 17, 20, 30, 40) = %.4f %.4f %.4f %.4f %.4f\n", kVname[x], fEps[x]->GetChisquare(), fEps[x]->GetNDF(),
               epsX17 / eps5 * fEps[x]->Eval(10) / fEps[0]->Eval(17), epsX17 / eps5 * fEps[x]->Eval(17) / fEps[0]->Eval(17), epsX17 / eps5 * fEps[x]->Eval(20) / fEps[0]->Eval(17),
               epsX17 / eps5 * fEps[x]->Eval(30) / fEps[0]->Eval(17), epsX17 / eps5 * fEps[x]->Eval(40) / fEps[0]->Eval(17));

    // output file
    TFile* fo = TFile::Open(Form("%s/x17_lines.root", kOutDir), "RECREATE");
    fMean.Write(); fRms.Write();
    for (int x = 0; x < kNVX; ++x) fEps[x]->Write();
    gMean.Write("gMean"); gRms.Write("gRms");
    for (int x = 0; x < kNVX; ++x) gEps[x].Write(Form("gEps_%s", kVname[x]));
    {
        TH1D hm("c6moments", "1 truncated mean, 2 truncated rms of the C6 line (17 +- 12 MeV)", 2, 0, 2);
        hm.SetBinContent(1, tm6); hm.SetBinContent(2, tr6); hm.Write();
    }
    c6.hee->Write("C6_ee"); c6.hm4->Write("C6_m4"); c5.hee->Write("C5V3_ee"); c5.hm4->Write("C5V3_m4");
    TTree* t = new TTree("lines", "X17 line per mass: slice of class-5 rows reweighted to X17 kinematics (S0, X2D M4 range)");
    double tm, teps[kNVX], tepse[kNVX], tR[kNVX], tmu, ts1, ts2, tf1, tmean, trms, tpg;
    long long tn;
    t->Branch("m", &tm, "m/D");
    t->Branch("eps", teps, Form("eps[%d]/D", kNVX));
    t->Branch("epsErr", tepse, Form("epsErr[%d]/D", kNVX));
    t->Branch("RX", tR, Form("RX[%d]/D", kNVX));
    t->Branch("dg_mu", &tmu, "dg_mu/D"); t->Branch("dg_s1", &ts1, "dg_s1/D"); t->Branch("dg_s2", &ts2, "dg_s2/D"); t->Branch("dg_f1", &tf1, "dg_f1/D");
    t->Branch("mean", &tmean, "mean/D"); t->Branch("rms", &trms, "rms/D"); t->Branch("pgen", &tpg, "pgen/D"); t->Branch("nrows", &tn, "nrows/L");
    for (int k = 0; k < kNM; ++k) {
        tm = kMlo + k * kMstep; tpg = pgen[k]; tn = sl[k][0].n;
        for (int x = 0; x < kNVX; ++x) {
            teps[x] = epsS[x][k]; tepse[x] = epsSe[x][k];
            // anchored on C6: eps_X(m) = eps_X(17, C6) x eps_slice(m) / eps_slice(17) (same variant)
            tR[x] = (epsS[0][k17] > 0) ? epsX17 * epsS[x][k] / epsS[0][k17] / eps5 : 0.;
            sl[k][x].hee->Write(); sl[k][x].hm4->Write();
        }
        tmu = dgs[k].mu; ts1 = dgs[k].s1; ts2 = dgs[k].s2; tf1 = dgs[k].f1; tmean = dgs[k].mean; trms = dgs[k].rms;
        t->Fill();
    }
    t->Write();
    TNamed info("info", Form("eps_X(17,C6) = %.6e; eps_5(C5 V3) = %.6e; R_X(17,C6) = %.5f; sigma_C6 = %.4f nb per unit BR; sigma_C5 = %.5f nb; "
                             "S0, X2D M4' range %.4f-%.4f, all M(ee); slices +-%.1f MeV; variants nom, DL, C2, NP",
                             epsX17, eps5, epsX17 / eps5, sigC6, sigC5, kM4lo, kM4hi, kHalf));
    info.Write();
    TH1D hc("consts", "1 epsX17, 2 eps5, 3 RX17, 4 sigC6, 5 sigC5, 6 epsSlice17 nom, 7 its error", 7, 0, 7);
    hc.SetBinContent(1, epsX17); hc.SetBinContent(2, eps5); hc.SetBinContent(3, epsX17 / eps5); hc.SetBinContent(4, sigC6);
    hc.SetBinContent(5, sigC5); hc.SetBinContent(6, epsS[0][k17]); hc.SetBinContent(7, epsSe[0][k17]);
    hc.Write();
    fo->Close();

    // table
    FILE* tb = fopen(Form("%s/x17ab_lines.txt", kTabDir), "w");
    fprintf(tb, "# x17ab_lines.txt -- X17 line model per mass (x17ab_line.C). Simulation only: eta set, owned rows (mode A), selection S0,\n");
    fprintf(tb, "# M4' = 0.99615 [m4_true + (m4 - m4_true)] in 0.5225-0.6100 GeV/c2, all M(ee).\n");
    fprintf(tb, "# C6 (m_X = 17.000 MeV/c2): sigma_vis = %.4f nb per unit BR(eta->pi+pi-X17) x BR(X17->ee); eps_X(17) = %.5e\n", sigC6, epsX17);
    fprintf(tb, "# C5 (V3): sigma_vis = %.5f nb; eps_5 = %.5e; R_X(17) = eps_X/eps_5 = %.4f\n", sigC5, eps5, epsX17 / eps5);
    fprintf(tb, "# validation at 17 MeV/c2 (C5 slice 16.5-17.5, reweighted, shifted) against C6:\n");
    fprintf(tb, "#   eps_slice(17)/eps_X(C6) = %.4f +- %.4f\n", epsS[0][k17] / epsX17, epsSe[0][k17] / epsX17);
    fprintf(tb, "#   line mean/rms C6 %.3f/%.3f, slice %.3f/%.3f MeV/c2; DG (mu, s1, s2, f1) C6 %.3f %.3f %.3f %.3f, slice %.3f %.3f %.3f %.3f\n",
            dg6.mean, dg6.rms, dgs[k17].mean, dgs[k17].rms, dg6.mu, dg6.s1, dg6.s2, dg6.f1, dgs[k17].mu, dgs[k17].s1, dgs[k17].s2, dgs[k17].f1);
    fprintf(tb, "#   shape chi2 C6 vs slice (1 MeV bins, 5-32 MeV/c2, non-empty bins) = %.1f / %d; M4' shape chi2 C6 vs slice %.1f / %d, C6 vs C5(V3) %.1f / %d\n",
            chiS, ndfS, chiM4b, ndfM2, chiM4a, ndfM);
    fprintf(tb, "# variants: nom = %s; DL = %s; C2 = %s; NP = %s\n", kVtext[0], kVtext[1], kVtext[2], kVtext[3]);
    fprintf(tb, "# R_X(m) = eps_X(17, C6) x eps_slice^var(m)/eps_slice^nom(17) / eps_5 (columns below: raw slices)\n");
    fprintf(tb, "# smooth (pol2 in m, 8-45 MeV/c2, used by the model): truncated mean - m = %.3f %+.4f m %+.6f m^2 (chi2 %.1f/%d); truncated rms = %.3f %+.4f m %+.6f m^2 (chi2 %.1f/%d)\n",
            fMean.GetParameter(0), fMean.GetParameter(1), fMean.GetParameter(2), fMean.GetChisquare(), fMean.GetNDF(), fRms.GetParameter(0),
            fRms.GetParameter(1), fRms.GetParameter(2), fRms.GetChisquare(), fRms.GetNDF());
    fprintf(tb, "# C6 truncated (17 +- 12 MeV/c2) mean %.3f rms %.3f MeV/c2\n", tm6, tr6);
    for (int x = 0; x < kNVX; ++x)
        fprintf(tb, "# smooth R_X %s (anchored on C6 at 17 through the NOMINAL slice efficiency: variant/nominal ratio kept at 17): m = 10 %.4f, 15 %.4f, 17 %.4f, 20 %.4f, 25 %.4f, 30 %.4f, 35 %.4f, 40 %.4f (fit chi2 %.1f/%d)\n", kVname[x],
                epsX17 / eps5 * fEps[x]->Eval(10) / fEps[0]->Eval(17), epsX17 / eps5 * fEps[x]->Eval(15) / fEps[0]->Eval(17), epsX17 / eps5 * fEps[x]->Eval(17) / fEps[0]->Eval(17),
                epsX17 / eps5 * fEps[x]->Eval(20) / fEps[0]->Eval(17), epsX17 / eps5 * fEps[x]->Eval(25) / fEps[0]->Eval(17),
                epsX17 / eps5 * fEps[x]->Eval(30) / fEps[0]->Eval(17), epsX17 / eps5 * fEps[x]->Eval(35) / fEps[0]->Eval(17),
                epsX17 / eps5 * fEps[x]->Eval(40) / fEps[0]->Eval(17), fEps[x]->GetChisquare(), fEps[x]->GetNDF());
    fprintf(tb, "# m[MeV]  rows  P_gen(slice)  eps_slice(nom)   err      R_X nom  R_X DL  R_X C2  R_X NP   mean   rms   DG: mu   s1    s2    f1\n");
    for (int k = 0; k < kNM; ++k) {
        double R[kNVX];
        for (int x = 0; x < kNVX; ++x) R[x] = (epsS[0][k17] > 0) ? epsX17 * epsS[x][k] / epsS[0][k17] / eps5 : 0.;
        fprintf(tb, "%6.1f %6lld %11.4e %12.5e %9.2e  %7.4f %7.4f %7.4f %7.4f  %6.2f %5.2f  %6.2f %5.2f %5.2f %5.3f\n", kMlo + k * kMstep, sl[k][0].n, pgen[k],
                epsS[0][k], epsSe[0][k], R[0], R[1], R[2], R[3], dgs[k].mean, dgs[k].rms, dgs[k].mu, dgs[k].s1, dgs[k].s2, dgs[k].f1);
    }
    fclose(tb);

    // figure 1: validation at 17 MeV
    {
        TCanvas c("c", "", 1400, 600);
        c.Divide(2, 1);
        c.cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1D* a = (TH1D*)c6.hee->Clone("a6"); TH1D* b = (TH1D*)sl[k17][0].hee->Clone("b17"); TH1D* b2 = (TH1D*)sl[k17][1].hee->Clone("b17dl");
        for (TH1D* h : {a, b, b2}) { h->SetDirectory(nullptr); h->Rebin(4); h->Scale(1. / h->Integral()); }
        a->SetLineColor(kBlack); a->SetMarkerStyle(20); a->SetMarkerColor(kBlack);
        b->SetLineColor(kRed + 1); b->SetLineWidth(2); b2->SetLineColor(kAzure + 2); b2->SetLineStyle(2); b2->SetLineWidth(2);
        a->GetXaxis()->SetRangeUser(0, 45);
        a->SetTitle(";M(e^{+}e^{-}) [MeV/c^{2}];fraction per 1 MeV/c^{2}");
        a->GetXaxis()->SetTitleSize(0.055); a->GetYaxis()->SetTitleSize(0.055); a->GetXaxis()->SetLabelSize(0.045); a->GetYaxis()->SetLabelSize(0.045);
        a->SetMaximum(1.55 * a->GetMaximum());
        a->Draw("E1"); b->Draw("HIST SAME"); b2->Draw("HIST SAME");
        TLegend lg(0.45, 0.70, 0.89, 0.89); lg.SetBorderSize(0); lg.SetTextSize(0.035); lg.SetFillStyle(0);
        lg.AddEntry(a, "X17 simulated (C6), m_{X} = 17.0", "lp");
        lg.AddEntry(b, "C5 slice 16.5-17.5, X17 kinematics", "l");
        lg.AddEntry(b2, "C5 slice, direct-like M(#pi#pi)", "l");
        lg.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.037);
        tx.DrawLatex(0.17, 0.85, "S0, M4' 0.5225-0.610");
        tx.DrawLatex(0.17, 0.80, Form("#chi^{2}/ndf = %.1f/%d", chiS, ndfS));
        tx.DrawLatex(0.17, 0.75, Form("#varepsilon_{slice}/#varepsilon_{C6} = %.3f #pm %.3f", epsS[0][k17] / epsX17, epsSe[0][k17] / epsX17));
        c.cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1D* m6 = (TH1D*)c6.hm4->Clone("m6"); TH1D* m5 = (TH1D*)c5.hm4->Clone("m5"); TH1D* ms = (TH1D*)sl[k17][0].hm4->Clone("ms");
        for (TH1D* h : {m6, m5, ms}) { h->SetDirectory(nullptr); h->Scale(1. / h->Integral()); }
        m6->SetLineColor(kBlack); m6->SetMarkerStyle(20); m5->SetLineColor(kBlue + 1); m5->SetLineWidth(2); ms->SetLineColor(kRed + 1); ms->SetLineWidth(2); ms->SetLineStyle(2);
        m6->SetTitle(";M'(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];fraction per 2.5 MeV/c^{2}");
        m6->GetXaxis()->SetTitleSize(0.055); m6->GetYaxis()->SetTitleSize(0.055); m6->GetXaxis()->SetLabelSize(0.045); m6->GetYaxis()->SetLabelSize(0.045);
        m6->GetXaxis()->SetNdivisions(505);
        m6->SetMaximum(1.3 * std::max(m6->GetMaximum(), m5->GetMaximum()));
        m6->Draw("E1"); m5->Draw("HIST SAME"); ms->Draw("HIST SAME");
        TLegend l2(0.50, 0.68, 0.89, 0.89); l2.SetBorderSize(0); l2.SetTextSize(0.037);
        l2.AddEntry(m6, "X17 (C6)", "lp"); l2.AddEntry(m5, "direct decay (C5, V3)", "l"); l2.AddEntry(ms, "C5 slice at 17, X17 kin.", "l"); l2.Draw();
        tx.DrawLatex(0.17, 0.85, Form("#chi^{2}/ndf C6 vs C5 = %.0f/%d", chiM4a, ndfM));
        c.SaveAs(Form("%s/x17ab_line_valid.pdf", kFigDir));
        c.SaveAs(Form("%s/x17ab_line_valid.png", kFigDir));
    }
    // figure 2: R_X(m) per variant and line widths
    {
        TCanvas c("c2", "", 1400, 600);
        c.Divide(2, 1);
        c.cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetGridy();
        TGraphErrors* g[kNVX];
        const int col[kNVX] = {kBlack, kAzure + 2, kRed + 1, kGreen + 2};
        double ymax = 0;
        for (int x = 0; x < kNVX; ++x) {
            g[x] = new TGraphErrors();
            for (int k = 0; k < kNM; ++k) {
                const double m = kMlo + k * kMstep;
                if (m < 8 || m > 45) continue;
                const double R = epsX17 * epsS[x][k] / epsS[0][k17] / eps5;
                const double Re = epsX17 * epsSe[x][k] / epsS[0][k17] / eps5;
                g[x]->SetPoint(g[x]->GetN(), m, R); g[x]->SetPointError(g[x]->GetN() - 1, 0, Re);
                ymax = std::max(ymax, R);
            }
            g[x]->SetLineColor(col[x]); g[x]->SetMarkerColor(col[x]); g[x]->SetMarkerStyle(20 + x); g[x]->SetMarkerSize(0.7);
        }
        g[0]->SetTitle(";m_{X} [MeV/c^{2}];R_{X} = #varepsilon_{X}/#varepsilon_{5}");
        g[0]->GetXaxis()->SetTitleSize(0.055); g[0]->GetYaxis()->SetTitleSize(0.055);
        g[0]->SetMinimum(0); g[0]->SetMaximum(1.3 * ymax);
        g[0]->Draw("APL");
        for (int x = 1; x < kNVX; ++x) g[x]->Draw("PL SAME");
        TGraph pc6; pc6.SetPoint(0, 17., epsX17 / eps5); pc6.SetMarkerStyle(29); pc6.SetMarkerSize(2.2); pc6.SetMarkerColor(kMagenta + 1); pc6.Draw("P SAME");
        TGraph gsm;
        for (int q = 0; q <= 74; ++q) { const double m = 8 + 0.5 * q; gsm.SetPoint(q, m, epsX17 / eps5 * fEps[0]->Eval(m) / fEps[0]->Eval(17)); }
        gsm.SetLineColor(kMagenta + 1); gsm.SetLineWidth(3); gsm.Draw("L SAME");
        TLegend lg(0.40, 0.16, 0.89, 0.42); lg.SetBorderSize(0); lg.SetTextSize(0.035);
        lg.AddEntry(g[0], "nominal (phase space, isotropic)", "lp"); lg.AddEntry(g[1], "direct-like M(#pi#pi)", "lp");
        lg.AddEntry(g[2], "e^{+} angle 1+cos^{2}", "lp"); lg.AddEntry(g[3], "pion angle sin^{2} kept", "lp"); lg.AddEntry(&pc6, "C6 (X17 sample)", "p"); lg.AddEntry(&gsm, "model: pol2, anchored on C6", "l");
        lg.Draw();
        c.cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetGridy();
        // truncated (m +- 12 MeV) mean - m and RMS of the slice lines with their pol2 curves (used by the model), C6 at 17
        TGraphErrors* gm = (TGraphErrors*)gMean.Clone("gm_draw");
        TGraphErrors* gr = (TGraphErrors*)gRms.Clone("gr_draw");
        gr->SetTitle(";m_{X} [MeV/c^{2}];line mean - m_{X}, line RMS [MeV/c^{2}]");
        gr->SetMinimum(-2); gr->SetMaximum(7);
        gr->SetMarkerStyle(20); gr->SetMarkerSize(0.7); gm->SetMarkerStyle(22); gm->SetMarkerSize(0.7); gm->SetMarkerColor(kRed + 1); gm->SetLineColor(kRed + 1);
        gr->GetXaxis()->SetTitleSize(0.055); gr->GetYaxis()->SetTitleSize(0.05);
        gr->Draw("AP"); gm->Draw("P SAME");
        TF1* fr2 = (TF1*)fRms.Clone("fr2"); TF1* fm2 = (TF1*)fMean.Clone("fm2");
        fr2->SetLineColor(kBlack); fm2->SetLineColor(kRed + 1); fr2->Draw("SAME"); fm2->Draw("SAME");
        TGraph* p6a = new TGraph(); p6a->SetPoint(0, 17, tr6); TGraph* p6c = new TGraph(); p6c->SetPoint(0, 17, tm6 - 17);
        for (TGraph* p : {p6a, p6c}) { p->SetMarkerStyle(29); p->SetMarkerSize(2.); p->SetMarkerColor(kMagenta + 1); p->Draw("P SAME"); }
        TLegend* l2 = new TLegend(0.40, 0.70, 0.89, 0.89); l2->SetBorderSize(0); l2->SetTextSize(0.035);
        l2->AddEntry(gr, "RMS (slices, #pm12 MeV/c^{2})", "p"); l2->AddEntry(gm, "mean - m_{X} (slices)", "p"); l2->AddEntry(p6a, "C6 (17 MeV/c^{2})", "p");
        l2->Draw();
        c.SaveAs(Form("%s/x17ab_line_eff.pdf", kFigDir));
        c.SaveAs(Form("%s/x17ab_line_eff.png", kFigDir));
    }
    printf("x17ab_line: done\n");
}
