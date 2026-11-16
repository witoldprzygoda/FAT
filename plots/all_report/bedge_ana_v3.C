// ============================================================================
// bedge_ana_v3.C  --  copy of bedge_ana.C for the ruling D38 (15 November 2026). The C2 gamma* factor of the per-slice
//                    estimator is the one the per-row QED weights of package c2recon give for its step (c2recon_bedge.C, density A):
//                    E1 0.868178, E2 0.853877 (all 0.859863) instead of the WP8 means 0.879 / 0.882; the visible yields use the
//                    window means of the same per-row weights (window B: E1 0.864249, E2 0.859664; all M4: E1 0.865095, E2 0.862860).
//                    Everything else is bedge_ana.C unchanged; inputs are read from the bedge directory, outputs go to bedge_v3/.
// bedge_ana.C  --  package bedge: f_eta(edge) as a measurement -- fits, errors, systematic variants
// ============================================================================
// Written 2026-11-14 (prefix bedge). Reads the arrays of bedge_fill.C and
//   1. reproduces WP6: f_eta(edge) = E_data / E_sim = 1141.9 / 2915.7 = 0.392 (fits of scan::FitEdgeStep, frozen
//      procedure of scan_bump.C), the WP6 MC error (100 replicas) and the WP6 18-fit model spread;
//   2. measures the statistical error with 200 event-group bootstrap replicas of the data (instead of the fit error x 1.2
//      of Part II) and the MC error with the 100 simulation replicas;
//   3. evaluates every systematic variant (families E1-E6, see the table kVars below) with the noise of each variant
//      difference measured on the same data replicas (bootstrap) and the same simulation replicas;
//   4. decomposes the simulated step into C2, C3 and the rest (the simulation fitted with the final data covariance,
//      i.e. the same linear estimator; exact decomposition), applies the HGeant gamma* angle correction of class C2
//      (WP8 factors, and their mean = recommended), and solves the two M(ee) classes for the bump eta scale a_D and the
//      conversion scale kappa (E_data,i = a (c_i C2_i + kappa C3_i + R_i), i = E1, E2);
//   5. per-period and per-fold values;
//   6. writes /hdd3/.../bedge/ana_results.txt (machine-readable) and bedge_figdata.json (for bedge_figs.py).
//
// Run (from /home/przygoda/HADES/FAT/FAT; about 2 GB, a few minutes):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_bedge", kTRUE)' \
//        'plots/all_report/bedge_ana.C+O'
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TMath.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "wp4_lib.h"
#include "bedge_common.h"

namespace bedgea {
using namespace bedge;
FILE* gL = nullptr;   // log
FILE* gR = nullptr;   // machine-readable results
#define PL(...) do { printf(__VA_ARGS__); if (gL) fprintf(gL, __VA_ARGS__); fflush(stdout); } while (0)

std::vector<double> DN, DR, SS, SV, SR, SW, SWV, SWR, HQ;
const double kBwD = (scan::kMhi - scan::kMlo) / scan::kNbin;   // = wp6 (FillMany and SimSpectra)
const double kRLo = scan::kMlo + kRM4lo * kBwD;                // sub-window low edge

// ---------------------------------------------------------------------------------------------------- slices
struct ESlice { double lo, hi, hr; bool useEdge = false; double edge = 0; };
inline int PB(double x) { return (int)std::lround((x - kPPlo) / kPPw); }
inline double HrRule(double lo) { return lo >= 0.39 - 1e-9 ? 0.035 : 0.045; }   // wp6 SliceCfg: 35 MeV from 0.39 on
std::vector<std::vector<ESlice>> gSets;
std::vector<std::string> gSetText;
std::vector<ESlice> gExtra;
inline ESlice MkS(double a, double b) { ESlice s; s.lo = a; s.hi = b; s.hr = HrRule(a); return s; }
void InitSets() {
    gSets.push_back({MkS(.30, .34), MkS(.34, .37), MkS(.37, .39), MkS(.39, .40)}); gSetText.push_back("nominal: 0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40");
    gSets.push_back({MkS(.30, .34), MkS(.34, .37), MkS(.37, .39), MkS(.39, .40), MkS(.40, .41)}); gSetText.push_back("+ 0.40-0.41");
    gSets.push_back({MkS(.30, .34), MkS(.34, .37), MkS(.37, .39)}); gSetText.push_back("without 0.39-0.40");
    gSets.push_back({MkS(.30, .32), MkS(.32, .34), MkS(.34, .355), MkS(.355, .37), MkS(.37, .38), MkS(.38, .39), MkS(.39, .395), MkS(.395, .40)});
    gSetText.push_back("8 half-width slices");
    gSets.push_back({MkS(.30, .37), MkS(.37, .40)}); gSetText.push_back("2 coarse slices 0.30-0.37, 0.37-0.40");
    gSets.push_back({MkS(.29, .34), MkS(.34, .37), MkS(.37, .39), MkS(.39, .40)}); gSetText.push_back("lower start 0.29");
    gSets.push_back({MkS(.31, .34), MkS(.34, .37), MkS(.37, .39), MkS(.39, .40)}); gSetText.push_back("lower start 0.31");
    gSets.push_back({MkS(.30, .33), MkS(.33, .36), MkS(.36, .38), MkS(.38, .40)}); gSetText.push_back("inner edges moved: 0.30-0.33-0.36-0.38-0.40");
    gSets.push_back({MkS(.30, .35), MkS(.35, .38), MkS(.38, .40)}); gSetText.push_back("3 slices 0.30-0.35-0.38-0.40");
    gExtra.push_back(MkS(.40, .41));
    ESlice c = MkS(.42, .46); c.useEdge = true; c.edge = 0.5150; c.hr = 0.045; gExtra.push_back(c);   // control slice (wp6 kCtrlEdge, half range 45 MeV)
}

// ---------------------------------------------------------------------------------------------------- spectra
void DCounts(int dc, int cls, int pa, int pb, double* u, double* p, double* m) {
    std::fill(u, u + kNM4, 0.); std::fill(p, p + kNM4, 0.); std::fill(m, m + kNM4, 0.);
    for (int e = 0; e < kNCls; ++e) {
        if (cls >= 0 && e != cls) continue;
        for (int q = pa; q < pb; ++q) {
            const double* a = &DN[IdxDN(dc, e, q, 0, 0)];
            const double* b = &DN[IdxDN(dc, e, q, 1, 0)];
            const double* c = &DN[IdxDN(dc, e, q, 2, 0)];
            for (int i = 0; i < kNM4; ++i) { u[i] += a[i]; p[i] += b[i]; m[i] += c[i]; }
        }
    }
}
void DCountsRep(int rc, int cls, int pa, int pb, int b, double* u, double* p, double* m) {
    std::fill(u, u + kNRM4, 0.); std::fill(p, p + kNRM4, 0.); std::fill(m, m + kNRM4, 0.);
    for (int e = 0; e < kNCls; ++e) {
        if (cls >= 0 && e != cls) continue;
        for (int q = pa; q < pb; ++q)
            for (int i = 0; i < kNRM4; ++i) {
                u[i] += DR[IdxDR(rc, e, q - kRPPlo, 0, i) + b];
                p[i] += DR[IdxDR(rc, e, q - kRPPlo, 1, i) + b];
                m[i] += DR[IdxDR(rc, e, q - kRPPlo, 2, i) + b];
            }
    }
}
scan::Spectra DSpec(int dc, int cls, int pa, int pb, const CBOpt& o, int clsRule) {
    std::vector<double> u(kNM4), p(kNM4), m(kNM4);
    DCounts(dc, -1, pa, pb, u.data(), p.data(), m.data());
    scan::Spectra all = MakeSpec(u.data(), p.data(), m.data(), kNM4, scan::kMlo, kBwD, o);
    if (cls < 0) return all;
    DCounts(dc, cls, pa, pb, u.data(), p.data(), m.data());
    scan::Spectra c = MakeSpec(u.data(), p.data(), m.data(), kNM4, scan::kMlo, kBwD, o);
    return clsRule == 0 ? SpreadCB(all, c) : c;
}
scan::Spectra DSpecRep(int rc, int cls, int pa, int pb, int b, const CBOpt& o, int clsRule) {
    std::vector<double> u(kNRM4), p(kNRM4), m(kNRM4);
    DCountsRep(rc, -1, pa, pb, b, u.data(), p.data(), m.data());
    scan::Spectra all = MakeSpec(u.data(), p.data(), m.data(), kNRM4, kRLo, kBwD, o);
    if (cls < 0) return all;
    DCountsRep(rc, cls, pa, pb, b, u.data(), p.data(), m.data());
    scan::Spectra c = MakeSpec(u.data(), p.data(), m.data(), kNRM4, kRLo, kBwD, o);
    return clsRule == 0 ? SpreadCB(all, c) : c;
}
const unsigned kMT12 = 0x1FF, kMC2 = 1u << gC2, kMC3 = 1u << gC3;
std::vector<double> SY(int sc, unsigned gm, int cls, int pa, int pb, double scale = 1.0) {
    std::vector<double> y(kNM4, 0.);
    for (int g = 0; g < kNG; ++g) {
        if (!(gm & (1u << g))) continue;
        for (int e = 0; e < kNCls; ++e) {
            if (cls >= 0 && e != cls) continue;
            for (int q = pa; q < pb; ++q) {
                const double* a = &SS[IdxSS(sc, g, e, q, 0)];
                for (int i = 0; i < kNM4; ++i) y[i] += scale * a[i];
            }
        }
    }
    return y;
}
std::vector<double> SYRep(int sc, int rg, int cls, int pa, int pb, int r, double scale = 1.0) {
    std::vector<double> y(kNRM4, 0.);
    for (int e = 0; e < kNCls; ++e) {
        if (cls >= 0 && e != cls) continue;
        for (int q = pa; q < pb; ++q)
            for (int i = 0; i < kNRM4; ++i) y[i] += scale * SR[IdxSR(sc, rg, e, q - kRPPlo, i) + r];
    }
    return y;
}

// ---------------------------------------------------------------------------------------------------- variants
struct Var {
    std::string tag, fam, text;
    int pair = -1;          // two variants with the same pair id: mean of |D+| and |D-|
    int counted = 1;        // 0: check only (listed, not counted)
    int r1 = 0;             // data-rejection test of WP7 (rule R1): 0 not applicable, 1 same bins, 2 other bins
    int dc = dS0, sc = sNom, set = 0, order = 2;
    double hrAll = -1, dLo = 0, dHi = 0;
    double sD = scan::kEdgeSigma, sS = scan::kEdgeSigma;
    int tmpl = 0;
    double scaleF = 1.0;
    int wgt = 0;            // 0: simulation fitted with its own Poisson-like weights (WP6), 1: with the final data covariance
    CBOpt cb;
    int clsRule = 0;        // 0: class CB by the spread rule (WP3/WP6), 1: the class's own 2 sqrt(N++ N--)
    double simScale = 1.0;
    int rc = 0;             // data replica configuration (0 S0, 1 pre, 2 S0best), -1: none
};
scan::StepCfg MakeCfg(const ESlice& s, const Var& v, bool sim) {
    scan::StepCfg c;
    c.mpLo = s.lo; c.mpHi = s.hi; c.order = v.order;
    c.halfRange = v.hrAll > 0 ? v.hrAll : s.hr;
    c.s = sim ? v.sS : v.sD;
    c.scale = scan::kScale * v.scaleF;
    if (s.useEdge) { c.useEdge = true; c.edge = s.edge * v.scaleF; }
    if (v.dLo != 0 || v.dHi != 0) {
        const double edge = s.useEdge ? c.edge : scan::EdgeMeasured(0.5 * (s.lo + s.hi), c.scale, c.scalePiPi);
        double lo = edge - c.halfRange - v.dLo, hi = edge + c.halfRange + v.dHi;
        int we;
        if (hi <= c.etaMean - 3.0 * c.etaSigma) we = 0;
        else if (hi < c.etaMean - c.etaSigma) { hi = c.etaMean - 3.0 * c.etaSigma; we = 0; }
        else { hi = std::max(hi, c.etaMean + 3.0 * c.etaSigma); we = 1; }
        c.lo = lo; c.hi = hi; c.withEta = we;
    }
    return c;
}
Tmpl MakeTmpl(int type, int cls, int pa, int pb) {
    Tmpl t; t.type = type;
    if (type == 1) {
        std::vector<double> h(kNHq, 0.);
        for (int e = 0; e < kNCls; ++e) {
            if (cls >= 0 && e != cls) continue;
            for (int q = pa; q < pb; ++q) for (int b = 0; b < kNHq; ++b) h[b] += HQ[IdxHq(e, q, b)];
        }
        for (int b = 0; b < kNHq; ++b) if (h[b] > 0) { t.q.push_back((b + 0.5) * 0.001); t.wq.push_back(h[b]); }
    }
    return t;
}

// ---------------------------------------------------------------------------------------------------- evaluation
struct SlR { double h = 0, eh = 0, chi2 = 0; int ndf = 0; double hs = 0, ehs = 0, chi2s = 0; int ndfs = 0; double lin[3] = {0, 0, 0}; bool ok = false; };
struct ClsR {
    double Ed = 0, vEd = 0, Es = 0, chi2 = 0; int ndf = 0; double lin[3] = {0, 0, 0}; bool ok = true;
    std::vector<SlR> sl; std::vector<std::vector<double>> Vd; std::vector<Tmpl> tD, tS; std::vector<scan::StepCfg> cD, cS;
};
struct Res { ClsR c[3]; bool ok() const { return c[0].ok && c[1].ok && c[2].ok; } };
struct Keep { std::vector<scan::StepFit> fd, fs; std::vector<scan::Spectra> d; std::vector<std::vector<double>> ys, y23; };

void EvaluateSlices(const Var& v, const std::vector<ESlice>& set, Res& R, Keep* keep = nullptr) {
    for (int k = 0; k < 3; ++k) {
        const int cls = k - 1;
        ClsR& C = R.c[k];
        C = ClsR();
        const int n = (int)set.size();
        C.sl.resize(n); C.Vd.resize(n); C.tD.resize(n); C.tS.resize(n); C.cD.resize(n); C.cS.resize(n);
        for (int i = 0; i < n; ++i) {
            const ESlice& s = set[i];
            const int pa = PB(s.lo), pb = PB(s.hi);
            C.cD[i] = MakeCfg(s, v, false); C.cS[i] = MakeCfg(s, v, true);
            C.tD[i] = MakeTmpl(v.tmpl, cls, pa, pb); C.tS[i] = C.tD[i];
            Tmpl* tD = v.tmpl ? &C.tD[i] : nullptr;
            Tmpl* tS = v.tmpl ? &C.tS[i] : nullptr;
            const scan::Spectra d = DSpec(v.dc, cls, pa, pb, v.cb, v.clsRule);
            const scan::StepFit fd = FitStep(d, C.cD[i], tD, nullptr, &C.Vd[i]);
            const std::vector<double> y = SY(v.sc, kMT12, cls, pa, pb, v.simScale);
            const scan::Spectra sp = SimSpec(y.data(), kNM4, scan::kMlo, kBwD);
            const scan::StepFit fs = v.wgt == 1 ? FitStep(sp, C.cS[i], tS, &C.Vd[i]) : FitStep(sp, C.cS[i], tS);
            SlR& S = C.sl[i];
            S.ok = fd.ok && fs.ok; S.h = fd.h; S.eh = fd.eh; S.chi2 = fd.chi2; S.ndf = fd.ndf; S.hs = fs.h; S.ehs = fs.eh; S.chi2s = fs.chi2; S.ndfs = fs.ndf;
            const unsigned gm[3] = {kMT12, kMC2, kMC3};
            for (int g = 0; g < 3; ++g) {
                const std::vector<double> yg = SY(v.sc, gm[g], cls, pa, pb, v.simScale);
                const scan::StepFit fl = FitStep(SimSpec(yg.data(), kNM4, scan::kMlo, kBwD), C.cS[i], tS, &C.Vd[i]);
                S.lin[g] = fl.h; C.lin[g] += fl.h;
                if (!fl.ok) S.ok = false;
            }
            C.Ed += fd.h; C.vEd += fd.eh * fd.eh; C.Es += fs.h; C.chi2 += fd.chi2; C.ndf += fd.ndf;
            if (!S.ok) C.ok = false;
            if (keep && k == 0) {
                keep->fd.push_back(fd); keep->fs.push_back(fs); keep->d.push_back(d); keep->ys.push_back(y);
                keep->y23.push_back(SY(v.sc, kMC2 | kMC3, cls, pa, pb, v.simScale));
            }
        }
    }
}
void Evaluate(const Var& v, Res& R, Keep* keep = nullptr) { EvaluateSlices(v, gSets[v.set], R, keep); }
// data replica b of configuration rc: E_data per class
void EvalBoot(const Var& v, Res& R, int b, double Ed[3]) {
    const std::vector<ESlice>& set = gSets[v.set];
    for (int k = 0; k < 3; ++k) {
        Ed[k] = 0;
        for (size_t i = 0; i < set.size(); ++i) {
            const int pa = PB(set[i].lo), pb = PB(set[i].hi);
            const scan::Spectra d = DSpecRep(v.rc, k - 1, pa, pb, b, v.cb, v.clsRule);
            const scan::StepFit fd = FitStep(d, R.c[k].cD[i], v.tmpl ? &R.c[k].tD[i] : nullptr);
            Ed[k] += fd.h;
        }
    }
}
// simulation replica r: E_sim (T12) per class and the linear decomposition (nominal data covariance)
void EvalRep(const Var& v, Res& R, int r, double Es[3], double lin[3][3]) {
    const std::vector<ESlice>& set = gSets[v.set];
    for (int k = 0; k < 3; ++k) {
        Es[k] = 0; lin[k][0] = lin[k][1] = lin[k][2] = 0;
        for (size_t i = 0; i < set.size(); ++i) {
            const int pa = PB(set[i].lo), pb = PB(set[i].hi);
            Tmpl* tS = v.tmpl ? &R.c[k].tS[i] : nullptr;
            for (int g = 0; g < 3; ++g) {
                const std::vector<double> y = SYRep(v.sc, g, k - 1, pa, pb, r, v.simScale);
                const scan::Spectra sp = SimSpec(y.data(), kNRM4, kRLo, kBwD);
                if (g == 0) {
                    const scan::StepFit fs = v.wgt == 1 ? FitStep(sp, R.c[k].cS[i], tS, &R.c[k].Vd[i]) : FitStep(sp, R.c[k].cS[i], tS);
                    Es[k] += fs.h;
                }
                const scan::StepFit fl = FitStep(sp, R.c[k].cS[i], tS, &R.c[k].Vd[i]);
                lin[k][g] += fl.h;
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------------- quantities
// C2 correction factors [set][class all, E1, E2]: 0 WP8 toy, 1 WP8 closure, 2 recommended = mean of the two
double CF[3][3] = {{0.880, 0.889, 0.868}, {0.891, 0.869, 0.896}, {0, 0, 0}};
const int kNQ = 14;
const char* const kQN[kNQ] = {"f_all", "f_E1", "f_E2", "fc_toy", "fc_clos", "fc_rec", "a_D", "kappa", "a_D_unc", "kappa_unc", "fc_E1", "fc_E2",
                              "a_k48", "a_k48_unc"};
const double kKapPeak = 0.48;   // conversion scale of the eta peak: conversion route / direct route = 0.237/0.497 (fix_numbers.txt)
struct Inp { double Ed[3], Es[3], lin[3][3]; };
void Quant(const Inp& x, double* q) {
    for (int k = 0; k < 3; ++k) q[k] = x.Ed[k] / x.Es[k];
    double sC2[3], sC3[3];
    for (int k = 0; k < 3; ++k) { sC2[k] = x.lin[k][1] / x.lin[k][0]; sC3[k] = x.lin[k][2] / x.lin[k][0]; }
    for (int j = 0; j < 3; ++j) {
        const double red = ((1 - CF[j][1]) * x.lin[1][1] + (1 - CF[j][2]) * x.lin[2][1]) / (x.lin[1][0] + x.lin[2][0]);
        q[3 + j] = x.Ed[0] / (x.Es[0] * (1 - red));
    }
    q[10] = x.Ed[1] / (x.Es[1] * (1 - (1 - CF[2][1]) * sC2[1]));
    q[11] = x.Ed[2] / (x.Es[2] * (1 - (1 - CF[2][2]) * sC2[2]));
    for (int corr = 0; corr < 2; ++corr) {
        const double c1 = corr == 0 ? CF[2][1] : 1.0, c2 = corr == 0 ? CF[2][2] : 1.0;
        const double A1 = x.Es[1] * (c1 * sC2[1] + 1 - sC2[1] - sC3[1]), B1 = x.Es[1] * sC3[1];
        const double A2 = x.Es[2] * (c2 * sC2[2] + 1 - sC2[2] - sC3[2]), B2 = x.Es[2] * sC3[2];
        const double rho = x.Ed[1] / x.Ed[2];
        const double kap = (rho * A2 - A1) / (B1 - rho * B2);
        const double a = x.Ed[2] / (A2 + kap * B2);
        q[6 + 2 * corr] = a; q[7 + 2 * corr] = kap;
    }
    // bump eta scale with the conversion part scaled by the peak's kappa (all M(ee)): E_data = a (c C2 + kappa C3 + R)
    for (int corr = 0; corr < 2; ++corr) {
        const double t = x.lin[1][0] + x.lin[2][0];
        const double c2 = corr == 0 ? ((CF[2][1] * x.lin[1][1] + CF[2][2] * x.lin[2][1]) / (x.lin[1][1] + x.lin[2][1])) : 1.0;
        const double s2 = (x.lin[1][1] + x.lin[2][1]) / t, s3 = (x.lin[1][2] + x.lin[2][2]) / t;
        q[12 + corr] = x.Ed[0] / (x.Es[0] * (c2 * s2 + kKapPeak * s3 + 1 - s2 - s3));
    }
}
Inp MakeInp(const Res& R) {
    Inp x;
    for (int k = 0; k < 3; ++k) { x.Ed[k] = R.c[k].Ed; x.Es[k] = R.c[k].Es; for (int g = 0; g < 3; ++g) x.lin[k][g] = R.c[k].lin[g]; }
    return x;
}
// statistical (fit) error of every quantity: numerical propagation of the independent E_data errors of the three classes
void StatFit(const Res& R, double* e) {
    const Inp x0 = MakeInp(R);
    double q0[kNQ]; Quant(x0, q0);
    double v[kNQ] = {0};
    for (int k = 0; k < 3; ++k) {
        Inp x = x0; const double d = std::sqrt(R.c[k].vEd);
        x.Ed[k] += d;
        double q[kNQ]; Quant(x, q);
        for (int j = 0; j < kNQ; ++j) v[j] += (q[j] - q0[j]) * (q[j] - q0[j]);
    }
    for (int j = 0; j < kNQ; ++j) e[j] = std::sqrt(v[j]);
}

// replica sets of one variant: q for every data replica (Es, lin nominal) and every simulation replica (Ed nominal)
struct RepQ { std::vector<std::array<double, kNQ>> b, r; };
void RunReplicas(const Var& v, Res& R, RepQ& Q, bool doBoot, bool doRep) {
    const Inp x0 = MakeInp(R);
    Q.b.clear(); Q.r.clear();
    if (doBoot && v.rc >= 0) {
        for (int b = 0; b < kNRepD; ++b) {
            Inp x = x0; EvalBoot(v, R, b, x.Ed);
            std::array<double, kNQ> q; Quant(x, q.data()); Q.b.push_back(q);
        }
    }
    if (doRep) {
        for (int r = 0; r < kNRepS; ++r) {
            Inp x = x0; EvalRep(v, R, r, x.Es, x.lin);
            std::array<double, kNQ> q; Quant(x, q.data()); Q.r.push_back(q);
        }
    }
}
double SD(const std::vector<double>& a) {
    if (a.size() < 2) return 0;
    double m = 0; for (double x : a) m += x; m /= a.size();
    double s = 0; for (double x : a) s += (x - m) * (x - m);
    return std::sqrt(s / (a.size() - 1));
}
double Mean(const std::vector<double>& a) { double m = 0; for (double x : a) m += x; return a.empty() ? 0 : m / a.size(); }
// half of the central 68 % interval (robust spread, WP7 "mc68")
double H68(std::vector<double> a) {
    if (a.size() < 5) return 0;
    std::sort(a.begin(), a.end());
    auto qt = [&](double p) { const double x = p * (a.size() - 1); const size_t i = (size_t)x; const double f = x - i; return i + 1 < a.size() ? a[i] * (1 - f) + a[i + 1] * f : a[i]; };
    return 0.5 * (qt(0.84134) - qt(0.15866));
}

}  // namespace bedgea

const char* const kOut3 = "/hdd3/przygoda/HADES/FAT_out/all_report/bedge_v3";
const double kCFW[2][2] = {{0.8642485658527286, 0.859664210021422}, {0.8650949589463703, 0.8628596663161261}};   // [window B, all][E1, E2]
void bedge_ana_v3(int quick = 0) {
    using namespace bedgea;
    TStopwatch total;
    gL = fopen(Form("%s/bedge_ana_v3.log", kOut3), "w");
    gR = fopen(Form("%s/ana_results_v3.txt", kOut3), "w");
    PL("# bedge_ana.C, 2026-11-14: f_eta(edge) = E_data/E_sim as a measurement\n");
    if (!ReadArr(Form("%s/data_N.bin", kOut), DN, SizeDN()) || !ReadArr(Form("%s/data_R.bin", kOut), DR, SizeDR()) ||
        !ReadArr(Form("%s/sim_S.bin", kOut), SS, SizeSS()) || !ReadArr(Form("%s/sim_V.bin", kOut), SV, SizeSS()) ||
        !ReadArr(Form("%s/sim_R.bin", kOut), SR, SizeSR()) || !ReadArr(Form("%s/sim_W.bin", kOut), SW, SizeSW()) ||
        !ReadArr(Form("%s/sim_WV.bin", kOut), SWV, SizeSW()) || !ReadArr(Form("%s/sim_WR.bin", kOut), SWR, SizeSWR()) ||
        !ReadArr(Form("%s/sim_hq.bin", kOut), HQ, SizeHq())) return;
    InitSets();
    CF[2][0] = 0.859863; CF[2][1] = 0.868178; CF[2][2] = 0.853877;   // D38: per-row QED weights through the bedge step (c2recon)
    fprintf(gR, "# ana_results_v3.txt of bedge_ana_v3.C (D38 per-row C2 factors); quantities: ");
    for (int j = 0; j < kNQ; ++j) fprintf(gR, "%d %s; ", j, kQN[j]);
    fprintf(gR, "\n# C2 factors (all, E1, E2): toy %.3f %.3f %.3f; closure %.3f %.3f %.3f; recommended (mean) %.4f %.4f %.4f\n", CF[0][0], CF[0][1],
            CF[0][2], CF[1][0], CF[1][1], CF[1][2], CF[2][0], CF[2][1], CF[2][2]);

    // ================================================================================================ 1. nominal and WP6 checks
    Var nom; nom.tag = "nominal"; nom.fam = "-"; nom.text = "WP6 frozen procedure";
    Res R0; Keep K0;
    Evaluate(nom, R0, &K0);
    PL("\n1. NOMINAL (S0, all M(ee)): E_data = %.4f +- %.4f, E_sim = %.4f, f = %.5f   [WP6: 1141.8584 +- 59.2243, 2915.6983, 0.39162]\n",
       R0.c[0].Ed, std::sqrt(R0.c[0].vEd), R0.c[0].Es, R0.c[0].Ed / R0.c[0].Es);
    PL("   E1: %.4f / %.4f = %.5f [WP6 636.0383 / 1734.7232];  E2: %.4f / %.4f = %.5f [WP6 505.5081 / 1171.8834]\n", R0.c[1].Ed, R0.c[1].Es,
       R0.c[1].Ed / R0.c[1].Es, R0.c[2].Ed, R0.c[2].Es, R0.c[2].Ed / R0.c[2].Es);
    // cross-check of the copied fitter against scan::FitEdgeStep
    {
        double maxd = 0;
        for (int i = 0; i < 4; ++i) {
            const scan::StepFit a = scan::FitEdgeStep(K0.d[i], R0.c[0].cD[i]);
            maxd = std::max(maxd, std::fabs(a.h - K0.fd[i].h));
        }
        PL("   copied fitter vs scan::FitEdgeStep (data, 4 slices): max |dh| = %.2e\n", maxd);
    }
    // sub-window reproduction (replica spectra are built on M4 0.40-0.62)
    {
        double Ed = 0; std::vector<double> u(kNRM4), p(kNRM4), m(kNRM4);
        for (int i = 0; i < 4; ++i) {
            const ESlice& s = gSets[0][i];
            std::vector<double> U(kNM4), Pp(kNM4), M(kNM4);
            DCounts(dS0, -1, PB(s.lo), PB(s.hi), U.data(), Pp.data(), M.data());
            const scan::Spectra d = MakeSpec(U.data() + kRM4lo, Pp.data() + kRM4lo, M.data() + kRM4lo, kNRM4, kRLo, kBwD);
            Ed += FitStep(d, R0.c[0].cD[i]).h;
        }
        PL("   the same fits on the sub-window 0.40-0.62: E_data = %.6f (full window %.6f)\n", Ed, R0.c[0].Ed);
    }
    // the simulated step split (data covariance: exact linear decomposition)
    for (int k = 0; k < 3; ++k)
        PL("   class %d: E_sim (WP6 weights) %.1f; with the data covariance: T12 %.1f = C2 %.1f + C3 %.1f + rest %.1f\n", k, R0.c[k].Es, R0.c[k].lin[0],
           R0.c[k].lin[1], R0.c[k].lin[2], R0.c[k].lin[0] - R0.c[k].lin[1] - R0.c[k].lin[2]);
    for (int k = 0; k < 3; ++k) {
        fprintf(gR, "NOM %d Ed %.6f eEd %.6f Es %.6f chi2 %.3f ndf %d linT12 %.6f linC2 %.6f linC3 %.6f\n", k, R0.c[k].Ed, std::sqrt(R0.c[k].vEd),
                R0.c[k].Es, R0.c[k].chi2, R0.c[k].ndf, R0.c[k].lin[0], R0.c[k].lin[1], R0.c[k].lin[2]);
        for (size_t i = 0; i < R0.c[k].sl.size(); ++i) {
            const SlR& S = R0.c[k].sl[i];
            fprintf(gR, "SLICE %d %zu %.3f %.3f h %.6f eh %.6f chi2 %.3f ndf %d hs %.6f ehs %.6f chi2s %.3f ndfs %d lin %.6f %.6f %.6f\n", k, i,
                    gSets[0][i].lo, gSets[0][i].hi, S.h, S.eh, S.chi2, S.ndf, S.hs, S.ehs, S.chi2s, S.ndfs, S.lin[0], S.lin[1], S.lin[2]);
        }
    }
    // extra slices (not in E): 0.40-0.41 and the control slice
    {
        Res Rx; EvaluateSlices(nom, gExtra, Rx);
        for (int k = 0; k < 3; ++k)
            for (size_t i = 0; i < gExtra.size(); ++i) {
                const SlR& S = Rx.c[k].sl[i];
                fprintf(gR, "EXTRA %d %zu %.3f %.3f h %.6f eh %.6f hs %.6f lin %.6f %.6f %.6f\n", k, i, gExtra[i].lo, gExtra[i].hi, S.h, S.eh, S.hs,
                        S.lin[0], S.lin[1], S.lin[2]);
                PL("   extra slice %.2f-%.2f class %d: data %.1f +- %.1f, sim %.1f (C2 %.1f, C3 %.1f)\n", gExtra[i].lo, gExtra[i].hi, k, S.h, S.eh, S.hs, S.lin[1], S.lin[2]);
            }
    }
    // WP6 18-fit joint model spread (order 2/3 x half range 35/45/55 x smearing 4.4/5.4/6.4; a varied slice fit with p <= 0.05
    // keeps its nominal h in data AND simulation)
    {
        double rmin = 1e9, rmax = -1e9; int kept = 0;
        const double f0 = R0.c[0].Ed / R0.c[0].Es;
        for (int ord = 2; ord <= 3; ++ord)
            for (int ih = 0; ih < 3; ++ih)
                for (int is = -1; is <= 1; ++is) {
                    Var v = nom; v.order = ord; v.hrAll = 0.035 + 0.010 * ih; v.sD = v.sS = scan::kEdgeSigma + 0.001 * is;
                    Res R; Evaluate(v, R);
                    double Ea = 0, Eb = 0; int nv = 0;
                    for (int s = 0; s < 4; ++s) {
                        const bool use = R.c[0].sl[s].ok && TMath::Prob(R.c[0].sl[s].chi2, R.c[0].sl[s].ndf) > 0.05;
                        Ea += use ? R.c[0].sl[s].h : R0.c[0].sl[s].h;
                        Eb += use ? R.c[0].sl[s].hs : R0.c[0].sl[s].hs;
                        if (use) ++nv;
                    }
                    if (nv == 0 || !(Eb > 0)) continue;
                    ++kept; rmin = std::min(rmin, Ea / Eb); rmax = std::max(rmax, Ea / Eb);
                }
        const double spread = 0.5 * (std::max(rmax, f0) - std::min(rmin, f0));
        PL("   WP6 18-fit model spread: %.5f over %d variations [WP6 0.02123]\n", spread, kept);
        fprintf(gR, "WP6GRID spread %.6f kept %d min %.6f max %.6f\n", spread, kept, rmin, rmax);
    }

    // ================================================================================================ 2. replicas of the nominal
    TStopwatch sw;
    RepQ Q0; RunReplicas(nom, R0, Q0, true, true);
    double q0[kNQ], e0[kNQ]; Quant(MakeInp(R0), q0); StatFit(R0, e0);
    PL("\n2. NOMINAL quantities (value, stat fit, stat bootstrap, bootstrap mean, MC, replica mean); %.0f s\n", sw.RealTime());
    for (int j = 0; j < kNQ; ++j) {
        std::vector<double> b, r;
        for (auto& x : Q0.b) b.push_back(x[j]);
        for (auto& x : Q0.r) r.push_back(x[j]);
        PL("   %-10s %9.5f  fit %8.5f  boot %8.5f [68 %8.5f] (mean %9.5f)  MC %8.5f [68 %8.5f] (mean %9.5f)\n", kQN[j], q0[j], e0[j], SD(b), H68(b), Mean(b), SD(r), H68(r), Mean(r));
        fprintf(gR, "NOMQ %s %.6f statfit %.6f statboot %.6f boot68 %.6f bootmean %.6f mc %.6f mc68 %.6f repmean %.6f\n", kQN[j], q0[j], e0[j], SD(b), H68(b), Mean(b), SD(r), H68(r), Mean(r));
    }
    // per-replica values of the main quantities (for correlations and the figures)
    {
        FILE* fq = fopen(Form("%s/nominal_replicas_v3.txt", kOut3), "w");
        fprintf(fq, "# data bootstrap replicas (b) and simulation replicas (r) of the nominal: Ed[3] Es[3] then the %d quantities\n", kNQ);
        for (size_t b = 0; b < Q0.b.size(); ++b) { fprintf(fq, "b %zu", b); for (int j = 0; j < kNQ; ++j) fprintf(fq, " %.6f", Q0.b[b][j]); fprintf(fq, "\n"); }
        for (size_t r = 0; r < Q0.r.size(); ++r) { fprintf(fq, "r %zu", r); for (int j = 0; j < kNQ; ++j) fprintf(fq, " %.6f", Q0.r[r][j]); fprintf(fq, "\n"); }
        fclose(fq);
    }
    // simulation replicas of E_sim alone (WP6 check: rms 287.79, mean 2877.3)
    {
        std::vector<double> es;
        for (int r = 0; r < kNRepS; ++r) { double Es[3], lin[3][3]; EvalRep(nom, R0, r, Es, lin); es.push_back(Es[0]); }
        PL("   E_sim replicas: rms %.2f, mean %.2f [WP6 287.79, 2877.3]\n", SD(es), Mean(es));
        fprintf(gR, "ESIMREP rms %.4f mean %.4f\n", SD(es), Mean(es));
    }

    if (quick) { PL("quick mode: stop after the nominal\n"); fclose(gR); fclose(gL); return; }
    // ================================================================================================ 3. variants
    std::vector<Var> vars;
    auto add = [&](Var v) { vars.push_back(v); };
    int pid = 0;
    {   // E1 edge model
        Var v = nom; v.fam = "E1";
        v.tag = "smear4.4"; v.text = "smearing 4.4 MeV (data and simulation)"; v.sD = v.sS = 0.0044; v.pair = pid; v.r1 = 1; add(v);
        v = nom; v.fam = "E1"; v.tag = "smear6.4"; v.text = "smearing 6.4 MeV"; v.sD = v.sS = 0.0064; v.pair = pid++; v.r1 = 1; add(v);
        v = nom; v.fam = "E1"; v.tag = "smearMatch"; v.text = "each sample its own peak width: data 5.4, simulation 6.1 MeV"; v.sS = 0.0061; v.r1 = 0; add(v);
        v = nom; v.fam = "E1"; v.tag = "tmplQ"; v.text = "step averaged over the simulated M(ee) of the bump classes (edge E(m, q))"; v.tmpl = 1; v.r1 = 1; add(v);
        v = nom; v.fam = "E1"; v.tag = "tmplBox"; v.text = "smeared box in M4^2 (phase-space shape) instead of the step"; v.tmpl = 2; v.r1 = 1; add(v);
        v = nom; v.fam = "E1"; v.tag = "pol1"; v.text = "background pol1"; v.order = 1; v.r1 = 1; add(v);
        v = nom; v.fam = "E1"; v.tag = "pol3"; v.text = "background pol3"; v.order = 3; v.r1 = 1; add(v);
        v = nom; v.fam = "E1"; v.tag = "hr35"; v.text = "fit range edge -+ 35 MeV in all slices"; v.hrAll = 0.035; v.pair = pid; v.r1 = 2; add(v);
        v = nom; v.fam = "E1"; v.tag = "hr55"; v.text = "fit range edge -+ 55 MeV in all slices"; v.hrAll = 0.055; v.pair = pid++; v.r1 = 2; add(v);
        v = nom; v.fam = "E1"; v.tag = "loExt"; v.text = "fit range extended by 10 MeV below"; v.dLo = 0.010; v.pair = pid; v.r1 = 2; add(v);
        v = nom; v.fam = "E1"; v.tag = "hiExt"; v.text = "fit range extended by 10 MeV above"; v.dHi = 0.010; v.pair = pid++; v.r1 = 2; add(v);
        v = nom; v.fam = "E1"; v.tag = "edgeUp"; v.text = "edge position x 1.0009 (eta-peak position +1 sigma)"; v.scaleF = 1.0009; v.pair = pid; v.r1 = 1; add(v);
        v = nom; v.fam = "E1"; v.tag = "edgeDn"; v.text = "edge position x 0.9991"; v.scaleF = 0.9991; v.pair = pid++; v.r1 = 1; add(v);
        v = nom; v.fam = "E1"; v.tag = "wgtData"; v.text = "simulation fitted with the final data covariance (same linear estimator)"; v.wgt = 1; add(v);
    }
    {   // E2 slices
        for (size_t s = 1; s < gSets.size(); ++s) {
            Var v = nom; v.fam = "E2"; v.set = (int)s; v.tag = Form("set%zu", s); v.text = gSetText[s]; add(v);
        }
    }
    {   // E3 mass transfer
        Var v = nom; v.fam = "E3";
        v.tag = "sUp"; v.text = kSC[sSUp].text; v.sc = sSUp; v.pair = pid; add(v);
        v.tag = "sDn"; v.text = kSC[sSDn].text; v.sc = sSDn; v.pair = pid++; add(v);
        v.tag = "sK0"; v.text = kSC[sSK0].text; v.sc = sSK0; v.pair = -1; add(v);
        v.tag = "r089"; v.text = kSC[sR089].text; v.sc = sR089; v.pair = pid; add(v);
        v.tag = "r110"; v.text = kSC[sR110].text; v.sc = sR110; v.pair = pid++; add(v);
        v.tag = "ppUp"; v.text = kSC[sPPUp].text; v.sc = sPPUp; v.pair = pid; add(v);
        v.tag = "ppDn"; v.text = kSC[sPPDn].text; v.sc = sPPDn; v.pair = pid++; add(v);
    }
    {   // E4 CB
        Var v = nom; v.fam = "E4";
        v.tag = "k0.95"; v.text = "CB x 0.95"; v.cb.k = 0.95; v.pair = pid; add(v);
        v.tag = "k1.05"; v.text = "CB x 1.05"; v.cb.k = 1.05; v.pair = pid++; add(v);
        v.tag = "k1.032"; v.text = "CB x 1.032 (data event mixing)"; v.cb.k = 1.032; v.pair = -1; add(v);
        v = nom; v.fam = "E4"; v.tag = "cb5"; v.text = "CB in 5 MeV cells"; v.cb.reb = 2; add(v);
        v = nom; v.fam = "E4"; v.tag = "cb20"; v.text = "CB in 20 MeV cells (no bootstrap: cells not aligned with the replica window)"; v.cb.reb = 8; v.rc = -1; add(v);
        v = nom; v.fam = "E4"; v.tag = "cbSum"; v.text = "CB = N++ + N-- (like-sign sum)"; v.cb.mode = 1; add(v);
        v = nom; v.fam = "E4"; v.tag = "clsOwnCB"; v.text = "class CB = its own 2 sqrt(N++ N--) (E1, E2 only)"; v.clsRule = 1; add(v);
    }
    {   // E5 physics weights / ownership
        Var v = nom; v.fam = "E5";
        v.tag = "wtsOff"; v.text = kSC[sWOff].text; v.sc = sWOff; add(v);
        v.tag = "bump1aY"; v.text = kSC[sB1aY].text; v.sc = sB1aY; add(v);
        v.tag = "bumpPS"; v.text = kSC[sBPS].text; v.sc = sBPS; add(v);
        v.tag = "modeB"; v.text = kSC[sModeB].text; v.sc = sModeB; add(v);
    }
    {   // E6 selection
        Var v = nom; v.fam = "E6"; v.tag = "pre"; v.text = "selection pre (data and simulation)"; v.dc = dPre; v.sc = sPre; v.rc = 1; add(v);
    }
    {   // checks (not counted)
        Var v = nom; v.fam = "CHK"; v.counted = 0;
        v.tag = "best"; v.text = "isBest == 1 in data and simulation"; v.dc = dBest; v.sc = sBest; v.rc = 2; add(v);
        v.tag = "bestData"; v.text = "isBest == 1 in data, simulation isBest >= 0 (WP7 convention)"; v.dc = dBest; v.sc = sNom; v.rc = 2; add(v);
    }
    PL("\n3. VARIANTS (%zu)\n", vars.size());
    // reference replica values per quantity
    for (const Var& v : vars) {
        TStopwatch s1;
        Res R; Evaluate(v, R);
        RepQ Q; RunReplicas(v, R, Q, true, true);
        double q[kNQ], e[kNQ]; Quant(MakeInp(R), q); StatFit(R, e);
        fprintf(gR, "VAR %s fam %s pair %d counted %d r1 %d chi2 %.4f ndf %d chi2_0 %.4f ndf_0 %d ok %d boot %d | %s\n", v.tag.c_str(), v.fam.c_str(), v.pair,
                v.counted, v.r1, R.c[0].chi2, R.c[0].ndf, R0.c[0].chi2, R0.c[0].ndf, R.ok() ? 1 : 0, Q.b.empty() ? 0 : 1, v.text.c_str());
        for (int k = 0; k < 3; ++k) {
            fprintf(gR, "VSL %s %d", v.tag.c_str(), k);
            for (size_t i = 0; i < R.c[k].sl.size(); ++i) fprintf(gR, " %.5f %.5f %.5f", R.c[k].sl[i].h, R.c[k].sl[i].eh, R.c[k].sl[i].hs);
            fprintf(gR, "\n");
        }
        PL("   %-10s %-3s f %.5f (E %.1f / %.1f) chi2 %.1f/%d", v.tag.c_str(), v.fam.c_str(), q[0], R.c[0].Ed, R.c[0].Es, R.c[0].chi2, R.c[0].ndf);
        for (int j = 0; j < kNQ; ++j) {
            std::vector<double> db, dr, b, r;
            for (size_t i = 0; i < Q.b.size(); ++i) { db.push_back(Q.b[i][j] - Q0.b[i][j]); b.push_back(Q.b[i][j]); }
            for (size_t i = 0; i < Q.r.size(); ++i) { dr.push_back(Q.r[i][j] - Q0.r[i][j]); r.push_back(Q.r[i][j]); }
            fprintf(gR, "VQ %s %s %.6f statfit %.6f statboot %.6f mc %.6f delta %.6f noiseboot %.6f noiserep %.6f nb68 %.6f nr68 %.6f\n", v.tag.c_str(), kQN[j], q[j], e[j],
                    SD(b), SD(r), q[j] - q0[j], SD(db), SD(dr), H68(db), H68(dr));
        }
        {
            std::vector<double> db, dr;
            for (size_t i = 0; i < Q.b.size(); ++i) db.push_back(Q.b[i][0] - Q0.b[i][0]);
            for (size_t i = 0; i < Q.r.size(); ++i) dr.push_back(Q.r[i][0] - Q0.r[i][0]);
            PL("  d %+.5f nb %.5f nr %.5f  (%.0f s)\n", q[0] - q0[0], SD(db), SD(dr), s1.RealTime());
        }
    }

    // ================================================================================================ 4. periods and folds
    PL("\n4. PERIODS (data period / simulation x L_p/L; fit stat errors) and FOLDS\n");
    {
        const double Lp[4] = {2809.453, 372.610, 1078.022, 1819.110};   // analysed L per period (WP2/WP7), sum 6079.2
        const double Csec[4] = {1.0241, 1.0, 1.0, 1.0};
        for (int p = 0; p < 4; ++p) {
            Var v = nom; v.dc = dPerI + p; v.simScale = Lp[p] / kL; v.rc = -1;
            Res R; Evaluate(v, R);
            double q[kNQ], e[kNQ]; Quant(MakeInp(R), q); StatFit(R, e);
            fprintf(gR, "PER %s Lp %.3f Csec %.4f", kDTag[dPerI + p], Lp[p], Csec[p]);
            for (int j = 0; j < kNQ; ++j) fprintf(gR, " %s %.6f %.6f", kQN[j], q[j], e[j]);
            fprintf(gR, "\n");
            PL("   %-7s f %.4f +- %.4f (fit); x C_sec %.4f; E1 %.4f E2 %.4f\n", kDTag[dPerI + p], q[0], e[0], q[0] * Csec[p], q[1], q[2]);
        }
        // folds: row fraction of the e+e- rows in M(pipi) 0.30-0.40
        double nS0 = 0, nF[2] = {0, 0};
        for (int e = 0; e < kNCls; ++e)
            for (int pp = PB(0.30); pp < PB(0.40); ++pp)
                for (int i = 0; i < kNM4; ++i) {
                    nS0 += DN[IdxDN(dS0, e, pp, 0, i)];
                    nF[0] += DN[IdxDN(dFold0, e, pp, 0, i)]; nF[1] += DN[IdxDN(dFold1, e, pp, 0, i)];
                }
        for (int fo = 0; fo < 2; ++fo) {
            Var v = nom; v.dc = dFold0 + fo; v.simScale = nF[fo] / nS0; v.rc = -1;
            Res R; Evaluate(v, R);
            double q[kNQ], e[kNQ]; Quant(MakeInp(R), q); StatFit(R, e);
            fprintf(gR, "FOLD %d frac %.6f", fo, nF[fo] / nS0);
            for (int j = 0; j < kNQ; ++j) fprintf(gR, " %s %.6f %.6f", kQN[j], q[j], e[j]);
            fprintf(gR, "\n");
            PL("   fold %d (row fraction %.4f): f %.4f +- %.4f\n", fo, nF[fo] / nS0, q[0], e[0]);
        }
    }

    // ================================================================================================ 5. window sums of the bump classes
    PL("\n5. SMASH (C2 + C3) window sums [raw rows] (nominal), MC from the replicas\n");
    for (int w = 0; w < kNWS; ++w)
        for (int k = 0; k < 3; ++k) {
            double s2 = 0, s3 = 0, v2 = 0, v3 = 0;
            std::vector<double> r2(kNRepS, 0.), r3(kNRepS, 0.);
            for (int e = 0; e < kNCls; ++e) {
                if (k > 0 && e != k - 1) continue;
                s2 += SW[IdxSW(sNom, gC2, e, w)]; s3 += SW[IdxSW(sNom, gC3, e, w)];
                v2 += SWV[IdxSW(sNom, gC2, e, w)]; v3 += SWV[IdxSW(sNom, gC3, e, w)];
                for (int r = 0; r < kNRepS; ++r) { r2[r] += SWR[IdxSWR(sNom, 0, e, w) + r]; r3[r] += SWR[IdxSWR(sNom, 1, e, w) + r]; }
            }
            std::vector<double> rs(kNRepS);
            for (int r = 0; r < kNRepS; ++r) rs[r] = r2[r] + r3[r];
            PL("   window %s class %d: C2 %.1f, C3 %.1f, sum %.1f (x 64/63 = %.1f), MC rms %.1f\n", w == 0 ? "B" : "all M4", k, s2, s3, s2 + s3, (s2 + s3) * kPT3, SD(rs));
            fprintf(gR, "WIN %s %d C2 %.4f C3 %.4f sum %.4f vC2 %.4f vC3 %.4f mcsum %.4f mcC2 %.4f mcC3 %.4f\n", w == 0 ? "B" : "all", k, s2, s3, s2 + s3, v2, v3, SD(rs),
                    SD(r2), SD(r3));
        }

    // window sums for every simulation configuration (for the yields of the variants)
    for (int c = 0; c < kNSCfg; ++c)
        for (int w = 0; w < kNWS; ++w) {
            fprintf(gR, "WINCFG %s %s", kSC[c].tag, w == 0 ? "B" : "all");
            for (int k = 0; k < 3; ++k) {
                double s2 = 0, s3 = 0;
                for (int e = 0; e < kNCls; ++e) { if (k > 0 && e != k - 1) continue; s2 += SW[IdxSW(c, gC2, e, w)]; s3 += SW[IdxSW(c, gC3, e, w)]; }
                fprintf(gR, " %.4f %.4f", s2, s3);
            }
            fprintf(gR, "\n");
        }
    // visible bump yields with MC errors from the SAME replicas as E_sim (the window sum and the step share simulated events)
    //   Y1 = f x (C2 + C3)_W;  Y2 = fc_rec x (c C2 + C3)_W;  Y3 = a_k48 x (c C2 + 0.48 C3)_W   [raw rows; x 64/63 in the note]
    {
        const Inp x0 = MakeInp(R0);
        for (int w = 0; w < kNWS; ++w)
            for (int k = 0; k < 3; ++k) {
                auto W = [&](int r, double& w2, double& w3) {
                    w2 = w3 = 0;
                    for (int e = 0; e < kNCls; ++e) {
                        if (k > 0 && e != k - 1) continue;
                        w2 += r < 0 ? SW[IdxSW(sNom, gC2, e, w)] : SWR[IdxSWR(sNom, 0, e, w) + r];
                        w3 += r < 0 ? SW[IdxSW(sNom, gC3, e, w)] : SWR[IdxSWR(sNom, 1, e, w) + r];
                    }
                };
                // class-mean C2 factor for the window
                double w2, w3; W(-1, w2, w3);
                double c2 = k == 0 ? 0 : kCFW[w][k - 1];   // D38: window mean of the per-row QED weights
                if (k == 0) {
                    double a2 = 0, b2 = 0;
                    for (int e = 0; e < 2; ++e) { a2 += kCFW[w][e] * SW[IdxSW(sNom, gC2, e, w)]; b2 += SW[IdxSW(sNom, gC2, e, w)]; }
                    c2 = a2 / b2;
                }
                auto Y = [&](const double* q, double v2, double v3, double* y) {
                    y[0] = q[k] * (v2 + v3);
                    y[1] = (k == 0 ? q[5] : (k == 1 ? q[10] : q[11])) * (c2 * v2 + v3);
                    y[2] = q[12] * (c2 * v2 + kKapPeak * v3);
                };
                double q[kNQ]; Quant(x0, q);
                double y0[3]; Y(q, w2, w3, y0);
                std::vector<double> yr[3], yb[3];
                for (int r = 0; r < kNRepS; ++r) {
                    Inp x = x0; EvalRep(nom, R0, r, x.Es, x.lin);
                    double qq[kNQ]; Quant(x, qq);
                    double v2, v3; W(r, v2, v3);
                    double y[3]; Y(qq, v2, v3, y);
                    for (int j = 0; j < 3; ++j) yr[j].push_back(y[j]);
                }
                for (size_t b = 0; b < Q0.b.size(); ++b) { double y[3]; Y(Q0.b[b].data(), w2, w3, y); for (int j = 0; j < 3; ++j) yb[j].push_back(y[j]); }
                fprintf(gR, "YIELD %s %d W2 %.4f W3 %.4f c2 %.5f", w == 0 ? "B" : "all", k, w2, w3, c2);
                for (int j = 0; j < 3; ++j) fprintf(gR, " Y%d %.4f stat %.4f mc %.4f mc68 %.4f", j + 1, y0[j], SD(yb[j]), SD(yr[j]), H68(yr[j]));
                fprintf(gR, "\n");
                PL("   yield %s class %d: Y1 %.0f +- %.0f +- %.0f (MC68 %.0f); Y2 %.0f +- %.0f (MC68); Y3 %.0f +- %.0f (MC68)  [raw rows]\n", w == 0 ? "B" : "all M4", k,
                   y0[0], SD(yb[0]), SD(yr[0]), H68(yr[0]), y0[1], H68(yr[1]), y0[2], H68(yr[2]));
            }
    }

    // ================================================================================================ 5b. background step projection
    // The step template is odd around the edge, so the odd (cubic) part of a smooth background that pol2 cannot absorb
    // projects onto h. Smoothed (wp4::SmoothFit, exp(pol3) on 0.40-0.70, per slice and class) and raw simulated G_D and G_C,
    // fitted with the final data covariance (linear): E(G), E(x G) with x = (M4 - 0.55)/0.1, so that any stage-A-type
    // background c_j (1 + b x)(G_D + kappa_B G_C) can be combined; and the raw rest of the eta classes (C4 + C5 + G_eta).
    PL("\n5b. BACKGROUND STEP PROJECTION (data covariance; rows per 2.5 MeV summed over the 4 slices)\n");
    {
        std::vector<double> edges; for (int i = 0; i <= 120; ++i) edges.push_back(0.40 + 0.0025 * i);
        const int i0 = (int)std::lround((0.40 - scan::kMlo) / 0.0025);
        for (int k = 0; k < 3; ++k) {
            double E[2][2][2] = {{{0}}};   // [GD, GC][smooth, raw][1, x]
            double Ealt[2][2] = {{0, 0}, {0, 0}};   // [GD, GC][exp(pol2), range 0.43-0.63]
            double Eoe = 0;
            fprintf(gR, "BKGSL %d", k);
            for (int i = 0; i < 4; ++i) {
                const ESlice& sl = gSets[0][i];
                const int pa = PB(sl.lo), pb = PB(sl.hi);
                double es[2][2][2] = {{{0}}};
                for (int t = 0; t < 2; ++t) {
                    const unsigned gm = t == 0 ? ((1u << gGDom) | (1u << gGDnom)) : ((1u << gGCom) | (1u << gGCnom));
                    const std::vector<double> y = SY(sNom, gm, k - 1, pa, pb);
                    std::vector<double> v(kNM4, 0.);
                    for (int g = 0; g < kNG; ++g) {
                        if (!(gm & (1u << g))) continue;
                        for (int e = 0; e < kNCls; ++e) { if (k > 0 && e != k - 1) continue; for (int q = pa; q < pb; ++q) for (int b = 0; b < kNM4; ++b) v[b] += SV[IdxSS(sNom, g, e, q, b)]; }
                    }
                    wp4::Vec cc(122, 0.), vv(122, 0.);
                    for (int b = 1; b <= 120; ++b) { cc[b] = y[i0 + b - 1]; vv[b] = v[i0 + b - 1]; }
                    const wp4::Smooth S = wp4::SmoothFit(cc, vv, edges);
                    // smoothing variants of the projection: exp(pol2) on 0.40-0.70; exp(pol3) on 0.43-0.63
                    {
                        const wp4::Smooth S2 = wp4::SmoothFit(cc, vv, edges, nullptr, 2);
                        std::vector<double> e2; for (int b = 12; b <= 92; ++b) e2.push_back(edges[b]);
                        wp4::Vec c3(83, 0.), v3(83, 0.);
                        for (int b = 1; b <= 80; ++b) { c3[b] = cc[b + 12]; v3[b] = vv[b + 12]; }
                        const wp4::Smooth S3 = wp4::SmoothFit(c3, v3, e2);
                        std::vector<double> y2(kNM4, 0.), y3(kNM4, 0.);
                        for (int b = 0; b < 120; ++b) { y2[i0 + b] = S2.Int(edges[b], edges[b + 1]); if (b >= 12 && b < 92) y3[i0 + b] = S3.Int(edges[b], edges[b + 1]); }
                        const double h2 = FitStep(SimSpec(y2.data(), kNM4, scan::kMlo, kBwD), R0.c[k].cS[i], nullptr, &R0.c[k].Vd[i]).h;
                        const double h3 = FitStep(SimSpec(y3.data(), kNM4, scan::kMlo, kBwD), R0.c[k].cS[i], nullptr, &R0.c[k].Vd[i]).h;
                        Ealt[t][0] += h2; Ealt[t][1] += h3;
                    }
                    std::vector<double> ys(kNM4, 0.), ysx(kNM4, 0.), yrx(kNM4, 0.);
                    for (int b = 0; b < kNM4; ++b) {
                        const double xc = scan::kMlo + (b + 0.5) * kBwD;
                        if (b >= i0 && b < i0 + 120) { double a0, a1; S.Int2(edges[b - i0], edges[b - i0 + 1], a0, a1); ys[b] = a0; ysx[b] = a1; }
                        yrx[b] = y[b] * TiltX(xc);
                    }
                    const std::vector<double>* src[2][2] = {{&ys, &ysx}, {&y, &yrx}};
                    for (int sr = 0; sr < 2; ++sr)
                        for (int xx = 0; xx < 2; ++xx) {
                            const double h = FitStep(SimSpec(src[sr][xx]->data(), kNM4, scan::kMlo, kBwD), R0.c[k].cS[i], nullptr, &R0.c[k].Vd[i]).h;
                            es[t][sr][xx] = h; E[t][sr][xx] += h;
                        }
                }
                const std::vector<double> yo = SY(sNom, (1u << gC4) | (1u << gC5) | (1u << gGeta), k - 1, pa, pb);
                const double ho = FitStep(SimSpec(yo.data(), kNM4, scan::kMlo, kBwD), R0.c[k].cS[i], nullptr, &R0.c[k].Vd[i]).h;
                Eoe += ho;
                fprintf(gR, " %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f", es[0][0][0], es[0][0][1], es[1][0][0], es[1][0][1], es[0][1][0], es[0][1][1], es[1][1][0], es[1][1][1], ho);
            }
            fprintf(gR, "\n");
            fprintf(gR, "BKG %d GDs %.4f GDsx %.4f GCs %.4f GCsx %.4f GDr %.4f GDrx %.4f GCr %.4f GCrx %.4f Oeta %.4f\n", k, E[0][0][0], E[0][0][1], E[1][0][0], E[1][0][1],
                    E[0][1][0], E[0][1][1], E[1][1][0], E[1][1][1], Eoe);
            PL("   class %d: G_D smooth %.1f (x-weighted %.1f), raw %.1f; G_C smooth %.1f (%.1f), raw %.1f; other eta (raw) %.1f; bump C2 + C3 %.1f\n", k,
               E[0][0][0], E[0][0][1], E[0][1][0], E[1][0][0], E[1][0][1], E[1][1][0], Eoe, R0.c[k].lin[1] + R0.c[k].lin[2]);
            PL("            smoothing variants: G_D exp(pol2) %.1f, range 0.43-0.63 %.1f; G_C exp(pol2) %.1f, range %.1f\n", Ealt[0][0], Ealt[0][1], Ealt[1][0], Ealt[1][1]);
            fprintf(gR, "BKGALT %d GDp2 %.4f GDr43 %.4f GCp2 %.4f GCr43 %.4f\n", k, Ealt[0][0], Ealt[0][1], Ealt[1][0], Ealt[1][1]);
            // MC noise of the raw rest (T12 - C2 - C3) projection, from the simulation replicas (nominal data covariance)
            {
                std::vector<double> rr;
                for (int r = 0; r < kNRepS; ++r) { double Es[3], lin[3][3]; EvalRep(nom, R0, r, Es, lin); rr.push_back(lin[k][0] - lin[k][1] - lin[k][2]); }
                PL("            raw rest (T12 - C2 - C3) %.1f, replica rms %.1f (MC noise of the background + other-eta projection)\n", R0.c[k].lin[0] - R0.c[k].lin[1] - R0.c[k].lin[2], SD(rr));
                fprintf(gR, "BKGNOISE %d rest %.4f rms %.4f\n", k, R0.c[k].lin[0] - R0.c[k].lin[1] - R0.c[k].lin[2], SD(rr));
            }
        }
    }

    // ================================================================================================ 6. figure data
    {
        FILE* fj = fopen(Form("%s/bedge_figdata_v3.json", kOut3), "w");
        fprintf(fj, "{\n\"f\": %.6f,\n\"slices\": [\n", q0[0]);
        for (size_t i = 0; i < K0.fd.size(); ++i) {
            const scan::StepFit& fd = K0.fd[i]; const scan::StepFit& fs = K0.fs[i];
            const scan::Spectra& d = K0.d[i];
            fprintf(fj, "{\"lo\": %.3f, \"hi\": %.3f, \"edge\": %.6f, \"edgeLo\": %.6f, \"edgeHi\": %.6f, \"h\": %.4f, \"eh\": %.4f, \"hs\": %.4f, \"chi2\": %.3f, \"ndf\": %d, \"chi2s\": %.3f,\n",
                    gSets[0][i].lo, gSets[0][i].hi, fd.edge, fd.edgeLo, fd.edgeHi, fd.h, fd.eh, fs.h, fd.chi2, fd.ndf, fs.chi2);
            auto arr = [&](const char* name, const std::vector<double>& a, bool last = false) {
                fprintf(fj, "\"%s\": [", name);
                for (size_t q = 0; q < a.size(); ++q) fprintf(fj, "%s%.6g", q ? ", " : "", a[q]);
                fprintf(fj, "]%s\n", last ? "" : ",");
            };
            arr("tx", fd.tx); arr("fy", fd.fy); arr("fy0", fd.fy0); arr("ty", fd.ty);
            arr("sfy", fs.fy); arr("sfy0", fs.fy0);
            std::vector<double> x, y, ey, cb, ys, y23;
            for (int b = 0; b < d.nbin; ++b) {
                const double xc = d.x(b);
                if (xc < 0.43 || xc > 0.60) continue;
                x.push_back(xc); y.push_back(d.sig[b]); ey.push_back(std::sqrt(d.sigVar[b])); cb.push_back(d.cb[b]);
                ys.push_back(K0.ys[i][b]); y23.push_back(K0.y23[i][b]);
            }
            arr("x", x); arr("y", y); arr("ey", ey); arr("cb", cb); arr("ys", ys); arr("y23", y23, true);
            fprintf(fj, "}%s\n", i + 1 < K0.fd.size() ? "," : "");
        }
        fprintf(fj, "]\n}\n");
        fclose(fj);
    }
    PL("\ndone in %.0f s\n", total.RealTime());
    fclose(gR); fclose(gL);
}
