// ============================================================================
// verify_res_numbers.C  --  independent re-derivation of the key numbers of the
//                           chapter "Results without simulation" (prefix res)
// ============================================================================
// Written by the verifier of the res chapter. It does NOT use scan_compact.root,
// ScanEngine.h or res_common.h: it reads the FAT outputs (tree scan_nt) directly,
// applies the selections with its own code, builds the M4 spectra and the
// CB = 2 sqrt(N++ N--) (10 MeV bins, linear interpolation to the 2.5 MeV bin
// centres, k = 1), and fits with ROOT TF1 / Minuit (Neyman chi2, errors
// sqrt(N + var CB), no CB covariance): the numbers must agree with the res
// chapter within a few per cent, the row counts exactly.
//   1. N and CB in the data peak window 0.53364-0.55416 and the P_main content S
//      (frozen Gaussian 0.5439 / 5.4 MeV + pol2, 0.5225-0.6100) for pre (isBest >= 0
//      and == 1), pre + OA > 2, peak nominal + vertex, DIRECT, PEAK, X17 region;
//      Y = S x 64/63 / 6079.2 nb^-1
//   2. ECAL-tagged eta T: own photon choice (closest to M(e+e-gamma) = 0.140 among
//      g1..g3), own m5c from four-vectors (photon energy scaled so that
//      M(e+e-gamma) = m_pi0), frozen 0.5463 / 6.2 MeV + pol2 in 0.49-0.62, CB in 20 MeV bins
//   3. side-band share of the peak content above 20 MeV/c2 (pre; six M(ee) classes,
//      CB per class, r from my own pol2 fit) and the X17-window D (pre + vertex)
//   4. edge step h in the M(pi+pi-) slice 0.37-0.39 (pre), own kinematic edge
//   5. pre rows per period (period mapping) and rows failing the FAT cuts (must be 0)
// Selections: pre = cor_ok == 1 && isBest >= 0 && eVertReco_z > -160 mm; peak nominal
//   + vertex = pre + oa_lab < 55 + M(pipi) > 0.30 + pip_isUsedVertex && pim_isUsedVertex;
//   DIRECT = pre + dz_pi < 40 mm, dz_pi = max(|pip_z - vz|, |pim_z - vz|); PEAK = DIRECT +
//   MM4 < 2.85; X17 region = pre + 15.5 <= M(ee) < 24.5 MeV/c2 + both pions in the vertex fit.
//   All variables with the energy-loss corrected momenta (*_cor columns).
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/output_pippim{epem,epep,emem}_exp.root
// Output : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/res_verify/verify_res_numbers.txt
// Run (from /home/przygoda/HADES/FAT/FAT, about 1 min, one process):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_verify_res", kTRUE)' 'plots/exp_report/verify_res_numbers.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"

namespace vres {

const char* const kDir = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530";
const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/res_verify";
const double kMeta = 0.547862, kMpi0 = 0.1349768;
const double kMeanP = 0.5439, kSigP = 0.0054;
const double kWlo = kMeanP - 1.9 * kSigP, kWhi = kMeanP + 1.9 * kSigP;     // 0.53364, 0.55416
const double kSBlo = 0.5225, kSBhi = 0.56530;
const double kL = 6079.2, kPT3 = 64.0 / 63.0;
const double kMeanT = 0.5463, kSigT = 0.0062;
const double kMeMeV = 0.51099895, kMpiMeV = 139.56995;
// M4 histograms
const int kNB = 240; const double kLO = 0.35, kHI = 0.95, kBW = 0.0025;
// m5c histograms
const int kNT = 160; const double kTLO = 0.40, kTHI = 0.80;

FILE* gF = nullptr;
void P(const char* fmt, ...) {
    va_list a, b; va_start(a, fmt); va_copy(b, a);
    vprintf(fmt, a); if (gF) vfprintf(gF, fmt, b);
    va_end(a); va_end(b);
}
double GF(double a, double b, double m, double s) {
    const double q = 1.0 / (std::sqrt(2.0) * s);
    return 0.5 * (TMath::Erf((b - m) * q) - TMath::Erf((a - m) * q));
}
// upper limit of M(pi+pi-e+e-) for eta -> pi+pi-pi0, pi0 -> e+e-gamma with the photon lost
// (m_ee -> 0): the photon takes its smallest energy in the eta frame, E = (E0 - p0)/2,
// M4^2 = m_eta^2 - 2 m_eta E
double M4max(double mpp) {
    const double E0 = (kMeta * kMeta + kMpi0 * kMpi0 - mpp * mpp) / (2 * kMeta);
    const double p0sq = E0 * E0 - kMpi0 * kMpi0;
    if (p0sq < 0) return -1;
    const double eg = 0.5 * (E0 - std::sqrt(p0sq));
    const double m2 = kMeta * kMeta - 2 * kMeta * eg;
    return m2 > 0 ? std::sqrt(m2) : -1;
}

struct Row {
    float m4, mee, mpipi, oa, oalab, mm4, vz, dzpi;
    int best, vtx, file, tagged; float m5c;
};
// selections (all on top of cor_ok == 1 && isBest >= 0)
const int kNS = 17;
const char* const kSN[kNS] = {"pre", "pre, isBest==1", "pre, OA>2", "peak nominal+vertex", "DIRECT", "PEAK", "X17 region",
                              "pre+vertex", "bump slice 0.37-0.39 (pre)", "pre M(ee) 0-5", "pre M(ee) 5-10", "pre M(ee) 10-20",
                              "pre M(ee) 20-50", "pre M(ee) 50-100", "pre M(ee) 100-270", "pre+vertex, isBest==1", "pre (tagged)"};
bool Pass(int s, const Row& r) {
    if (!(r.vz > -160.f)) return false;
    switch (s) {
        case 0: return true;
        case 1: return r.best == 1;
        case 2: return r.oa > 2.f;
        case 3: return r.oalab < 55.f && r.mpipi > 0.30f && r.vtx;
        case 4: return r.dzpi < 40.f;
        case 5: return r.dzpi < 40.f && r.mm4 < 2.85f;
        case 6: return r.mee >= 0.0155f && r.mee < 0.0245f && r.vtx;
        case 7: return r.vtx;
        case 8: return r.mpipi >= 0.37f && r.mpipi < 0.39f;
        case 9: return r.mee < 0.005f;
        case 10: return r.mee >= 0.005f && r.mee < 0.010f;
        case 11: return r.mee >= 0.010f && r.mee < 0.020f;
        case 12: return r.mee >= 0.020f && r.mee < 0.050f;
        case 13: return r.mee >= 0.050f && r.mee < 0.100f;
        case 14: return r.mee >= 0.100f && r.mee < 0.270f;
        case 15: return r.vtx && r.best == 1;
        case 16: return r.tagged;
    }
    return false;
}

// CB = 2 sqrt(N++ N--) in wide bins of reb fine bins, linearly interpolated to the fine-bin centres
struct CB { std::vector<double> c, v, W, Wv; int reb = 4, nw = 0; };
CB MakeCB(const TH1D* pp, const TH1D* mm, int reb) {
    CB o; o.reb = reb;
    const int nb = pp->GetNbinsX(); o.nw = nb / reb;
    o.W.assign(o.nw, 0); o.Wv.assign(o.nw, 0);
    for (int j = 0; j < o.nw; ++j) {
        double a = 0, b = 0;
        for (int q = 0; q < reb; ++q) { a += pp->GetBinContent(j * reb + q + 1); b += mm->GetBinContent(j * reb + q + 1); }
        o.W[j] = 2 * std::sqrt(a * b); o.Wv[j] = a + b;
    }
    o.c.assign(nb, 0); o.v.assign(nb, 0);
    for (int i = 0; i < nb; ++i) {
        const double t = (i + 0.5) / reb - 0.5;   // position in units of wide bins, relative to the first wide-bin centre
        int j0, j1; double f;
        if (t <= 0) { j0 = j1 = 0; f = 0; }
        else if (t >= o.nw - 1) { j0 = j1 = o.nw - 1; f = 0; }
        else { j0 = (int)t; j1 = j0 + 1; f = t - j0; }
        o.c[i] = ((1 - f) * o.W[j0] + f * o.W[j1]) / reb;
        o.v[i] = ((1 - f) * (1 - f) * o.Wv[j0] + f * f * o.Wv[j1]) / (reb * reb);
    }
    return o;
}
// integral of the interpolated CB over [a, b) (fractional fine bins) and its variance (wide bins independent)
void CBInt(const CB& cb, double lo, double bw, double a, double b, double& val, double& var, std::vector<double>* wj = nullptr) {
    std::vector<double> w(cb.nw, 0.);
    val = 0;
    for (int i = 0; i < (int)cb.c.size(); ++i) {
        const double x1 = lo + i * bw, x2 = x1 + bw;
        const double u = std::max(0., std::min(x2, b) - std::max(x1, a)) / bw;
        if (u <= 0) continue;
        val += u * cb.c[i];
        const double t = (i + 0.5) / cb.reb - 0.5;
        int j0, j1; double f;
        if (t <= 0) { j0 = j1 = 0; f = 0; } else if (t >= cb.nw - 1) { j0 = j1 = cb.nw - 1; f = 0; } else { j0 = (int)t; j1 = j0 + 1; f = t - j0; }
        w[j0] += u * (1 - f) / cb.reb; w[j1] += u * f / cb.reb;
    }
    var = 0;
    for (int j = 0; j < cb.nw; ++j) var += w[j] * w[j] * cb.Wv[j];
    if (wj) *wj = w;
}

struct Fit { double S = 0, eS = 0, chi2 = 0; int ndf = 0; double c[3] = {0, 0, 0}, x0 = 0, dx = 1; double h = 0, eh = 0; bool ok = false;
    double Pol(double x) const { const double t = (x - x0) / dx; return c[0] + c[1] * t + c[2] * t * t; } };

// Gaussian (frozen) + pol2 chi2 fit of y +- e in [lo, hi) (bin centres)
Fit FitPeak(const std::vector<double>& y, const std::vector<double>& e, double hlo, double bw, double mean, double sig, double lo, double hi) {
    const int nb = (int)y.size();
    TH1D h("hfit", "", nb, hlo, hlo + nb * bw); h.SetDirectory(nullptr);
    for (int i = 0; i < nb; ++i) { h.SetBinContent(i + 1, y[i]); h.SetBinError(i + 1, e[i] > 0 ? e[i] : 1.0); }
    Fit r; r.x0 = 0.5 * (lo + hi); r.dx = 0.5 * (hi - lo);
    const double x0 = r.x0, dx = r.dx;
    TF1 f("fpk", [=](double* x, double* p) {
        const double t = (x[0] - x0) / dx;
        return p[0] * GF(x[0] - 0.5 * bw, x[0] + 0.5 * bw, mean, sig) + p[1] + p[2] * t + p[3] * t * t; }, lo, hi - 1e-7, 4);
    f.SetParameters(1000, h.GetBinContent(h.FindBin(x0)), 0, 0);
    TFitResultPtr fr = h.Fit(&f, "RQS0N");
    if (!fr.Get() || !fr->IsValid()) return r;
    r.S = fr->Parameter(0); r.eS = fr->ParError(0); r.chi2 = fr->Chi2(); r.ndf = fr->Ndf();
    for (int k = 0; k < 3; ++k) r.c[k] = fr->Parameter(k + 1);
    r.ok = true;
    return r;
}
// edge step in one M(pi+pi-) slice: pol2 + eta Gaussian + h x averaged smeared step
Fit FitStep(const std::vector<double>& y, const std::vector<double>& e, double mpLo, double mpHi, double lo, double hi) {
    const int nb = (int)y.size();
    TH1D h("hst", "", nb, kLO, kHI); h.SetDirectory(nullptr);
    for (int i = 0; i < nb; ++i) { h.SetBinContent(i + 1, y[i]); h.SetBinError(i + 1, e[i] > 0 ? e[i] : 1.0); }
    std::vector<double> edges;
    for (int k = 0; k < 21; ++k) edges.push_back((kMeanP / kMeta) * M4max((mpLo + (k + 0.5) * (mpHi - mpLo) / 21) / 0.9970));
    Fit r; r.x0 = 0.5 * (lo + hi); r.dx = 0.5 * (hi - lo);
    const double x0 = r.x0, dx = r.dx;
    TF1 f("fst", [=](double* x, double* p) {
        const double t = (x[0] - x0) / dx;
        double st = 0;
        for (double ed : edges) st += 0.5 * TMath::Erfc((x[0] - ed) / (std::sqrt(2.0) * kSigP));
        st /= edges.size();
        return p[0] + p[1] * t + p[2] * t * t + p[3] * GF(x[0] - 0.5 * kBW, x[0] + 0.5 * kBW, kMeanP, kSigP) + p[4] * st; }, lo, hi - 1e-7, 5);
    f.SetParameters(h.GetBinContent(h.FindBin(x0)), 0, 0, 100, 100);
    TFitResultPtr fr = h.Fit(&f, "RQS0N");
    if (!fr.Get() || !fr->IsValid()) return r;
    r.h = fr->Parameter(4); r.eh = fr->ParError(4); r.chi2 = fr->Chi2(); r.ndf = fr->Ndf(); r.ok = true;
    return r;
}

} // namespace vres

void verify_res_numbers() {
    using namespace vres;
    TStopwatch sw;
    gSystem->mkdir(kOut, kTRUE);
    gF = fopen(Form("%s/verify_res_numbers.txt", kOut), "w");
    P("# verify_res_numbers.txt -- independent re-derivation (FAT scan_nt -> own selections, CB, TF1 fits)\n");

    TH1D* h[kNS][3];
    TH1D* ht[3];
    for (int s = 0; s < kNS; ++s) for (int k = 0; k < 3; ++k) { h[s][k] = new TH1D(Form("h%d_%d", s, k), "", kNB, kLO, kHI); h[s][k]->SetDirectory(nullptr); }
    for (int k = 0; k < 3; ++k) { ht[k] = new TH1D(Form("ht_%d", k), "", kNT, kTLO, kTHI); ht[k]->SetDirectory(nullptr); }
    // exact window counts and per-group accumulators (unlike-sign) for the side-band subtraction
    double nW[kNS][3] = {}, nSB[kNS][3] = {};
    std::vector<std::pair<int, int>> grp[kNS];   // (rows in W, rows in SB) per event group, unlike-sign
    double nWt[3] = {};
    std::vector<double> fscW;   // photon scale factor of tagged unlike-sign pre rows in the m5c eta window
    long nPer[4] = {0, 0, 0, 0}, badTrig = 0, badVz = 0, nRows[3] = {0, 0, 0}, nTagRows = 0;
    long nGroups = 0;

    const char* fn[3] = {"output_pippimepem_exp.root", "output_pippimepep_exp.root", "output_pippimemem_exp.root"};
    for (int k = 0; k < 3; ++k) {
        TFile f(Form("%s/%s", kDir, fn[k]), "READ");
        TTree* t = (TTree*)f.Get("scan_nt");
        t->SetBranchStatus("*", 0);
        const char* bn[] = {"cor_ok", "isBest", "eVertReco_z", "trigbit", "file_id", "runnumber", "event", "m_pippimepem_cor", "m_ee_cor",
                            "m_pippim_cor", "oa", "oa_pippim_epem_lab_cor", "mm_pippimepem_cor", "pip_isUsedVertex", "pim_isUsedVertex",
                            "pip_z", "pim_z", "ep_p_cor", "ep_theta", "ep_phi", "em_p_cor", "em_theta", "em_phi", "pip_p_cor", "pip_theta",
                            "pip_phi", "pim_p_cor", "pim_theta", "pim_phi", "g1_e", "g1_theta", "g1_phi", "g2_e", "g2_theta", "g2_phi",
                            "g3_e", "g3_theta", "g3_phi"};
        const int nbn = sizeof(bn) / sizeof(bn[0]);
        std::vector<float> v(nbn, 0.f);
        for (int b = 0; b < nbn; ++b) { t->SetBranchStatus(bn[b], 1); t->SetBranchAddress(bn[b], &v[b]); }
        enum { COR, BEST, VZ, TRIG, FID, RUN, EVT, M4, MEE, MPP, OA, OAL, MM4, PIPV, PIMV, PIPZ, PIMZ,
               EPP, EPT, EPF, EMP, EMT, EMF, PPP, PPT, PPF, PMP, PMT, PMF, G1E, G1T, G1F, G2E, G2T, G2F, G3E, G3T, G3F };
        const Long64_t n = t->GetEntries();
        float lf = -1, lr = -1, le = -1;
        std::vector<int> aW(kNS, 0), aS(kNS, 0);
        auto flush = [&]() {
            for (int s = 0; s < kNS; ++s) { if (aW[s] || aS[s]) grp[s].push_back(std::make_pair(aW[s], aS[s])); aW[s] = aS[s] = 0; }
        };
        const double d2r = TMath::DegToRad();
        auto p4 = [&](double p, double th, double ph, double m) {
            TLorentzVector q; q.SetXYZM(p * std::sin(th * d2r) * std::cos(ph * d2r), p * std::sin(th * d2r) * std::sin(ph * d2r), p * std::cos(th * d2r), m);
            return q; };
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            nRows[k]++;
            if ((int)v[TRIG] != 8192) badTrig++;
            if (!(v[VZ] > -500.f)) badVz++;
            if (k == 0 && (v[FID] != lf || v[RUN] != lr || v[EVT] != le)) { flush(); lf = v[FID]; lr = v[RUN]; le = v[EVT]; nGroups++; }
            if (v[COR] != 1.f || !(v[BEST] >= 0.f)) continue;
            Row r;
            r.m4 = v[M4]; r.mee = v[MEE]; r.mpipi = v[MPP]; r.oa = v[OA]; r.oalab = v[OAL]; r.mm4 = v[MM4]; r.vz = v[VZ];
            r.dzpi = std::max(std::fabs(v[PIPZ] - v[VZ]), std::fabs(v[PIMZ] - v[VZ]));
            r.best = (int)std::lround(v[BEST]); r.vtx = (v[PIPV] == 1.f && v[PIMV] == 1.f); r.file = (int)std::lround(v[FID]);
            // ECAL tag: good photon (g_e > 0) with M(e+e-gamma) closest to 0.140
            r.tagged = 0; r.m5c = -1;
            {
                const TLorentzVector ee = p4(v[EPP], v[EPT], v[EPF], kMeMeV) + p4(v[EMP], v[EMT], v[EMF], kMeMeV);
                const TLorentzVector pp = p4(v[PPP], v[PPT], v[PPF], kMpiMeV) + p4(v[PMP], v[PMT], v[PMF], kMpiMeV);
                double best = 1e9; TLorentzVector g; bool have = false; double meeg = -1;
                for (int j = 0; j < 3; ++j) {
                    const float e = v[G1E + 3 * j];
                    if (!(e > 0)) continue;
                    const TLorentzVector q = p4(e, v[G1T + 3 * j], v[G1F + 3 * j], 0.0);
                    const double m = (ee + q).M() / 1000.0;
                    if (std::fabs(m - 0.140) < best) { best = std::fabs(m - 0.140); g = q; have = true; meeg = m; }
                }
                if (have && r.mee < 0.10f && meeg >= 0.119 && meeg < 0.173) {
                    // scale the photon so that M(e+e-gamma) = m_pi0: (ee + f g)^2 = m_pi0^2
                    const double mpi0 = 1000.0 * kMpi0;
                    const double fsc = (mpi0 * mpi0 - ee.M2()) / (2.0 * ee.Dot(g));
                    r.tagged = 1; r.m5c = (ee + pp + fsc * g).M() / 1000.0;
                    if (k == 0 && r.vz > -160.f && r.m5c >= kMeanT - 1.9 * kSigT && r.m5c < kMeanT + 1.9 * kSigT) fscW.push_back(fsc);
                }
            }
            if (k == 0 && r.vz > -160.f && r.m4 == r.m4) {
                const int d = r.file;
                const int p = (d >= 34 && d <= 46) ? 0 : (d <= 48 ? 1 : (d <= 53 ? 2 : (d >= 59 && d <= 67 ? 3 : -1)));
                if (p >= 0) nPer[p]++;
                if (r.tagged) nTagRows++;
            }
            if (!(r.m4 == r.m4)) continue;
            const bool inW = r.m4 >= kWlo && r.m4 < kWhi;
            const bool inS = (r.m4 >= kSBlo && r.m4 < kWlo) || (r.m4 >= kWhi && r.m4 < kSBhi);
            for (int s = 0; s < kNS; ++s) {
                if (!Pass(s, r)) continue;
                if (s == 16) {   // tagged eta in m5c
                    ht[k]->Fill(r.m5c);
                    if (r.m5c >= kMeanT - 1.9 * kSigT && r.m5c < kMeanT + 1.9 * kSigT) nWt[k] += 1;
                    continue;
                }
                h[s][k]->Fill(r.m4);
                if (inW) { nW[s][k] += 1; if (k == 0) aW[s]++; }
                if (inS) { nSB[s][k] += 1; if (k == 0) aS[s]++; }
            }
        }
        if (k == 0) flush();
        f.Close();
    }
    P("rows read: us %ld, pp %ld, mm %ld; unlike-sign event groups %ld; rows with trigbit != 8192: %ld; rows with vz <= -500: %ld\n",
      nRows[0], nRows[1], nRows[2], nGroups, badTrig, badVz);
    P("pre unlike-sign rows (cor_ok, finite M4) per period: I %ld, IIA %ld, IIB %ld, III %ld  [res: 2049702, 291000, 830907, 1251729]\n",
      nPer[0], nPer[1], nPer[2], nPer[3]);

    // ---- 1. peak fits ----
    struct Ref { int s; double N, CB, S, eS, Y; };
    const Ref ref[] = {{0, 77513, 7949, 4671, 367, 0.7805}, {1, 62759, 6308, 4199, 330, 0.7017}, {2, 60935, 7793, 3632, 328, 0.6070},
                       {3, 45433, 2510, 3949, 276, 0.6598}, {4, 62969, 4784, 4492, 327, 0.7507}, {5, 56532, 4066, 4364, 310, 0.7293},
                       {6, 5258, 102, 440, 92, 0.0735}, {7, 55947, 3279, 4388, 307, 0.7332}};
    P("\n1. ETA PEAK (own CB and TF1 fit; res values in brackets)\n");
    P("%-22s | %7s %7s | %7s %6s | %6s %5s %8s | %8s | %s\n", "selection", "N_W", "[res]", "CB_W", "[res]", "S", "stat", "chi2/ndf", "Y [nb]", "[res S +- stat, Y]");
    std::vector<Fit> pk(kNS);
    std::vector<CB> cbs(kNS);
    for (int s = 0; s < 16; ++s) cbs[s] = MakeCB(h[s][1], h[s][2], 4);
    for (const Ref& q : ref) {
        const int s = q.s;
        std::vector<double> y(kNB), e(kNB);
        for (int i = 0; i < kNB; ++i) { const double nn = h[s][0]->GetBinContent(i + 1); y[i] = nn - cbs[s].c[i]; e[i] = std::sqrt(std::max(nn, 1.0) + cbs[s].v[i]); }
        pk[s] = FitPeak(y, e, kLO, kBW, kMeanP, kSigP, 0.5225, 0.6100);
        double cbw, cbv; CBInt(cbs[s], kLO, kBW, kWlo, kWhi, cbw, cbv);
        const double Y = pk[s].S * kPT3 / kL;
        P("%-22s | %7.0f %7.0f | %7.1f %6.0f | %6.0f %5.0f %4.1f/%-3d | %8.4f | [%4.0f +- %3.0f, %.4f]  dS/S = %+.3f\n", kSN[s], nW[s][0], q.N, cbw, q.CB,
          pk[s].S, pk[s].eS, pk[s].chi2, pk[s].ndf, Y, q.S, q.eS, q.Y, (pk[s].S - q.S) / q.S);
    }
    P("   Y check with the res S: pre 4671 x 64/63 / 6079.2 = %.4f nb\n", 4671 * kPT3 / kL);
    for (int s : {0, 5}) {   // CB x 1.03 and x 0.97
        double sk[2];
        for (int q = 0; q < 2; ++q) {
            const double kf = q == 0 ? 1.03 : 0.97;
            std::vector<double> y(kNB), e(kNB);
            for (int i = 0; i < kNB; ++i) { const double nn = h[s][0]->GetBinContent(i + 1); y[i] = nn - kf * cbs[s].c[i]; e[i] = std::sqrt(std::max(nn, 1.0) + kf * kf * cbs[s].v[i]); }
            sk[q] = FitPeak(y, e, kLO, kBW, kMeanP, kSigP, 0.5225, 0.6100).S - pk[s].S;
        }
        P("   %-10s S(CB x 1.03) - S = %+.1f, S(CB x 0.97) - S = %+.1f  [res: pre +1/-1, PEAK -1/+1]\n", kSN[s], sk[0], sk[1]);
    }

    // ---- 2. tagged eta ----
    {
        const CB c = MakeCB(ht[1], ht[2], 8);
        std::vector<double> y(kNT), e(kNT);
        for (int i = 0; i < kNT; ++i) { const double nn = ht[0]->GetBinContent(i + 1); y[i] = nn - c.c[i]; e[i] = std::sqrt(std::max(nn, 1.0) + c.v[i]); }
        const double bw = (kTHI - kTLO) / kNT;
        const Fit f = FitPeak(y, e, kTLO, bw, kMeanT, kSigT, 0.49, 0.62);
        double cbw, cbv; CBInt(c, kTLO, bw, kMeanT - 1.9 * kSigT, kMeanT + 1.9 * kSigT, cbw, cbv);
        P("\n2. ECAL-TAGGED ETA (pre, isBest >= 0): tagged unlike-sign rows (any m5c) %ld [res 683381]; N in the m5c window %.0f [res 7420], CB %.1f [res 253]\n",
          nTagRows, nWt[0], cbw);
        P("   T = %.0f +- %.0f, chi2/ndf %.1f/%d [res 3314 +- 105, 37.7/48]; Y_T = %.4f nb [res 0.5537]\n", f.S, f.eS, f.chi2, f.ndf, f.S * kPT3 / kL);
        P("   S/T (pre, isBest >= 0) with my S and T: %.3f [res 1.409]\n", pk[0].S / f.S);
        std::vector<double> q = fscW; std::sort(q.begin(), q.end());
        double mean = 0; for (double x : q) mean += x; mean /= std::max<size_t>(1, q.size());
        P("   photon energy scale factor f (E_gamma x f gives M(e+e-gamma) = m_pi0), rows in the m5c eta window: median %.3f, mean %.3f (n = %zu)\n",
          q.empty() ? 0. : q[q.size() / 2], mean, q.size());
    }

    // ---- 3. side-band subtraction ----
    auto ratioWSB = [&](const Fit& f) {
        double w = 0, b = 0; const int n = 2000;
        for (int q = 0; q < n; ++q) { const double x = kSBlo + (q + 0.5) * (kSBhi - kSBlo) / n; const double val = f.Pol(x);
            if (x >= kWlo && x < kWhi) w += val; else b += val; }
        return w / b; };
    const double fW = GF(kWlo, kWhi, kMeanP, kSigP), fSB = GF(kSBlo, kWlo, kMeanP, kSigP) + GF(kWhi, kSBhi, kMeanP, kSigP);
    auto sideband = [&](int s, double r, double& D, double& eD) {
        double cw, cwv, cs, csv; std::vector<double> ww, ws;
        CBInt(cbs[s], kLO, kBW, kWlo, kWhi, cw, cwv, &ww);
        double c1, c1v, c2, c2v; std::vector<double> w1, w2;
        CBInt(cbs[s], kLO, kBW, kSBlo, kWlo, c1, c1v, &w1); CBInt(cbs[s], kLO, kBW, kWhi, kSBhi, c2, c2v, &w2);
        cs = c1 + c2;
        double vcb = 0;
        for (int j = 0; j < cbs[s].nw; ++j) { const double cc = ww[j] - r * (w1[j] + w2[j]); vcb += cc * cc * cbs[s].Wv[j]; }
        double vg = 0; for (const auto& g : grp[s]) { const double a = g.first - r * g.second; vg += a * a; }
        D = (nW[s][0] - cw) - r * (nSB[s][0] - cs); eD = std::sqrt(vg + vcb);
        (void)csv;
    };
    {
        const double r = ratioWSB(pk[0]);
        P("\n3. SIDE BAND (pre): own r = %.4f [res 0.9163]; fW %.4f fSB %.4f\n", r, fW, fSB);
        double a = 0, a2 = 0, c = 0, c2 = 0;
        for (int s = 9; s <= 14; ++s) {
            double D, eD; sideband(s, r, D, eD);
            const double S = D / (fW - r * fSB), eS = eD / (fW - r * fSB);
            double Dr, eDr; sideband(s, 0.9163, Dr, eDr);
            P("   %-18s N_W %6.0f  S = %6.0f +- %4.0f   (with the res r = 0.9163: %6.0f +- %4.0f)\n", kSN[s], nW[s][0], S, eS,
              Dr / (fW - 0.9163 * fSB), eDr / (fW - 0.9163 * fSB));
            if (s >= 12) { a += S; a2 += eS * eS; } else { c += S; c2 += eS * eS; }
        }
        const double t = a + c, fsh = a / t, efsh = std::sqrt(c * c * a2 + a * a * c2) / (t * t);
        P("   share above 20 MeV/c2: %.3f +- %.3f  [res 0.234 +- 0.043]; sum of the classes %.0f [res 4798]\n", fsh, efsh, t);
        // X17 window
        std::vector<double> y(kNB), e(kNB);
        for (int i = 0; i < kNB; ++i) { const double nn = h[7][0]->GetBinContent(i + 1); y[i] = nn - cbs[7].c[i]; e[i] = std::sqrt(std::max(nn, 1.0) + cbs[7].v[i]); }
        const Fit fb = FitPeak(y, e, kLO, kBW, kMeanP, kSigP, 0.5225, 0.6100);
        const double rx = ratioWSB(fb);
        double D, eD; sideband(6, rx, D, eD);
        double cw, cwv; CBInt(cbs[6], kLO, kBW, kWlo, kWhi, cw, cwv);
        P("   X17 window (pre + vertex, 15.5-24.5 MeV/c2): own r = %.4f [res 0.9148]; N_W %.0f [res 5258], CB_W %.1f [res 101.7], D = %.1f +- %.1f [res 481.5 +- 99.0]\n",
          rx, nW[6][0], cw, D, eD);
        P("   like-sign rows in W for the X17 window: N++ %.0f, N-- %.0f [res 14, 223]\n", nW[6][1], nW[6][2]);
    }

    // ---- 4. edge step, slice 0.37-0.39, pre ----
    {
        std::vector<double> y(kNB), e(kNB);
        for (int i = 0; i < kNB; ++i) { const double nn = h[8][0]->GetBinContent(i + 1); y[i] = nn - cbs[8].c[i]; e[i] = std::sqrt(std::max(nn, 1.0) + cbs[8].v[i]); }
        const Fit f = FitStep(y, e, 0.37, 0.39, 0.4600, 0.5600);
        P("\n4. EDGE STEP (pre, M(pipi) 0.37-0.39, range 0.4600-0.5600, pol2 + eta + step): h = %.1f +- %.1f, chi2/ndf %.1f/%d [res 354.8 +- 28.0, 32.5/35]\n",
          f.h, f.eh, f.chi2, f.ndf);
        P("   own predicted edge at the slice centre: %.4f [res 0.5048]\n", (kMeanP / kMeta) * M4max(0.38 / 0.9970));
    }
    P("\nrun time %.0f s\n", sw.RealTime());
    fclose(gF); gF = nullptr;
}
