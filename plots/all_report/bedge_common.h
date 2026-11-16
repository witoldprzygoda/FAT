// ============================================================================
// bedge_common.h  --  shared layer of the package "bedge" (the eta bump measured through its kinematic edge)
// ============================================================================
// Written 2026-11-14 (prefix bedge). Used by bedge_fill.C (row pass), bedge_ana.C (fits, systematics, conversions)
// and bedge_toys.C (pseudo-data closure).
//
// Contents
//   1. binning and the configuration lists (data and simulation)
//   2. binary array I/O (the row pass writes flat double arrays; layout in <out>/layout.txt)
//   3. spectra: scan::Spectra (ScanEngine.h) built from summed counts, with the CB rule of ScanEngine::FillMany
//      (2 sqrt(N++ N--) in 10 MeV cells, interpolated), the class spread rule of wp6_fill.C (SpreadCB, copied) and
//      the variants k x CB, 5 / 20 MeV cells, like-sign sum
//   4. the edge-step fit: a copy of scan::FitEdgeStep (ScanEngine.h, frozen procedure of scan_bump.C) with three
//      additions that do not change the nominal result: (a) a custom step template, (b) a FIXED covariance matrix
//      (the simulation fitted with the final data covariance = the same linear estimator), (c) the final covariance
//      is returned. Nominal reproduction checked against scan::FitEdgeStep in bedge_ana.C.
// ============================================================================
#ifndef BEDGE_COMMON_H
#define BEDGE_COMMON_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "defs.h"

namespace bedge {
using namespace alldefs;

const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/bedge";

// ---------------------------------------------------------------------------------------------------------------
// 1. binning
// ---------------------------------------------------------------------------------------------------------------
// M4: 240 x 2.5 MeV, 0.35-0.95 (= ScanEngine kNbin, kMlo, kMhi). Data bin index as FillMany ((m - lo)/bw with
// bw = (hi - lo)/nbin), simulation as wp6_fill.C ((m - kMlo)/kBinW).
const int kNM4 = 240;
// M(pipi) on the data scale: 40 x 5 MeV from 0.27 (edges rounded to 1e-3 so that 0.30, 0.34, ... are the literal
// doubles of wp6_fill.C; the data compare float values as wp6_fill.C does)
const int kNPP = 40;
const double kPPlo = 0.27, kPPw = 0.005;
inline double PPEdge(int k) { return std::round((kPPlo + kPPw * k) * 1000.) / 1000.; }
// M(ee) classes: 0 E1 (< 10 MeV), 1 E2 (>= 10 MeV), 2 no value (NaN; only in "all")
const int kNCls = 3;
// replica sub-window: M4 bins [20, 108) = 0.40-0.62 (88 bins, aligned with the 10 MeV CB cells of the full grid),
// M(pipi) bins [2, 40) = 0.28-0.47
const int kRM4lo = 20, kNRM4 = 88;
const int kRPPlo = 2, kNRPP = 38;
const int kNRepD = 200;   // data bootstrap (Poisson(1) per event group, own seed)
const int kNRepS = 100;   // simulation replicas (WP0 generator and seed: replica r = hists_v1 replica r)
const uint64_t kSeedSim = 0x5157A11ULL;
const uint64_t kSeedData = 0xBED6E0DA7AULL;

// ---- counter-based Poisson(1) replica weights: copy of wp0_fill.C / wp6_fill.C -----------------------------------
inline uint64_t Mix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
inline int Pois1(uint64_t h) {
    const double u = (double)(h >> 11) * (1.0 / 9007199254740992.0);
    double p = 0.36787944117144233, c = p;
    int k = 0;
    while (u >= c && k < 30) { ++k; p /= k; c += p; }
    return k;
}
inline void GroupPois(uint64_t seed, int tree, int grp, int n, double* k) {
    const uint64_t base = Mix(seed ^ Mix(((uint64_t)(uint32_t)tree << 32) | (uint64_t)(uint32_t)grp));
    for (int r = 0; r < n; ++r) k[r] = Pois1(Mix(base + 0x632BE59BD9B4E019ULL * (uint64_t)(r + 1)));
}

// ---- data configurations -------------------------------------------------------------------------------------
enum DCfg { dS0 = 0, dPre, dBest, dPerI, dPerIIA, dPerIIB, dPerIII, dFold0, dFold1, kNDCfg };
const char* const kDTag[kNDCfg] = {"S0", "pre", "S0best", "perI", "perIIA", "perIIB", "perIII", "fold0", "fold1"};
const int kNDRep = 3;   // replicas for S0, pre, S0best

// ---- simulation configurations -------------------------------------------------------------------------------
struct SCfgDef { const char* tag; int sel; int best; int phys; double s, r, spp; int mode; const char* text; };
enum SCfg { sNom = 0, sSUp, sSDn, sSK0, sR089, sR110, sF, sWOff, sB1aY, sBPS, sPre, sBest, sPPUp, sPPDn, sModeB, kNSCfg };
const SCfgDef kSC[kNSCfg] = {
    {"nom", kSelS0, 0, kPhysNom, kS, kR, kSPiPi, 0, "nominal: S0, mode A, w x V3 (108), s = 0.99615, r = 1, M(pipi) x 0.9952"},
    {"sUp", kSelS0, 0, kPhysNom, kS + kSerr, kR, kSPiPi, 0, "s = 0.99705"},
    {"sDn", kSelS0, 0, kPhysNom, kS - kSerr, kR, kSPiPi, 0, "s = 0.99525"},
    {"sK0", kSelS0, 0, kPhysNom, kSK0, kR, kSPiPi, 0, "s = 0.99521 (K0S-based)"},
    {"r089", kSelS0, 0, kPhysNom, kS, kRlo, kSPiPi, 0, "r = 0.89 (data-like resolution)"},
    {"r110", kSelS0, 0, kPhysNom, kS, kRhi, kSPiPi, 0, "r = 1.10"},
    {"f", kSelS0, 0, kPhysNom, kS + kSerr, kRlo, kSPiPi, 0, "s = 0.99705 and r = 0.89 (WP4 scenario f)"},
    {"wtsOff", kSelS0, 0, kPhysOff, kS, kR, kSPiPi, 0, "physics weights off (0)"},
    {"bump1aY", kSelS0, 0, kPhysBump1aY, kS, kR, kSPiPi, 0, "eta -> 3pi Dalitz plot 1 + aY (116)"},
    {"bumpPS", kSelS0, 0, kPhysBumpPS, kS, kR, kSPiPi, 0, "eta -> 3pi phase space (100)"},
    {"pre", kSelPre, 0, kPhysNom, kS, kR, kSPiPi, 0, "selection pre"},
    {"best", kSelS0, 1, kPhysNom, kS, kR, kSPiPi, 0, "S0 and isBest == 1"},
    {"ppUp", kSelS0, 0, kPhysNom, kS, kR, kSPiPi + 0.0007, 0, "M(pipi) x 0.9959"},
    {"ppDn", kSelS0, 0, kPhysNom, kS, kR, kSPiPi - 0.0007, 0, "M(pipi) x 0.9945"},
    {"modeB", kSelS0, 0, kPhysNom, kS, kR, kSPiPi, 1, "ownership mode B"}};
// simulation groups of the spectra (the T12 sum is their sum)
enum SGrp { gC2 = 0, gC3, gC4, gC5, gGeta, gGDom, gGDnom, gGCom, gGCnom, kNG };
const char* const kGTag[kNG] = {"C2", "C3", "C4", "C5", "Geta", "GDom", "GDnom", "GCom", "GCnom"};
// replica groups
enum RGrp { rT12 = 0, rC2, rC3, kNRG };
// window sums: 0 = B (0.4550-0.5225), 1 = all M4; all M(pipi)
const int kNWS = 2;
const int kNHq = 150;   // M(ee) 0-0.150 in 1 MeV, C2 + C3 near the edge (nominal configuration)

// ---- flat indices ------------------------------------------------------------------------------------------------
inline size_t IdxDN(int c, int e, int p, int k, int i) { return ((((size_t)c * kNCls + e) * kNPP + p) * 3 + k) * kNM4 + i; }
inline size_t SizeDN() { return (size_t)kNDCfg * kNCls * kNPP * 3 * kNM4; }
inline size_t IdxDR(int c, int e, int p, int k, int i) { return (((((size_t)c * kNCls + e) * kNRPP + p) * 3 + k) * kNRM4 + i) * kNRepD; }
inline size_t SizeDR() { return (size_t)kNDRep * kNCls * kNRPP * 3 * kNRM4 * kNRepD; }
inline size_t IdxSS(int c, int g, int e, int p, int i) { return ((((size_t)c * kNG + g) * kNCls + e) * kNPP + p) * kNM4 + i; }
inline size_t SizeSS() { return (size_t)kNSCfg * kNG * kNCls * kNPP * kNM4; }
inline size_t IdxSR(int c, int g, int e, int p, int i) { return (((((size_t)c * kNRG + g) * kNCls + e) * kNRPP + p) * kNRM4 + i) * kNRepS; }
inline size_t SizeSR() { return (size_t)kNSCfg * kNRG * kNCls * kNRPP * kNRM4 * kNRepS; }
inline size_t IdxSW(int c, int g, int e, int w) { return (((size_t)c * kNG + g) * kNCls + e) * kNWS + w; }
inline size_t SizeSW() { return (size_t)kNSCfg * kNG * kNCls * kNWS; }
inline size_t IdxSWR(int c, int g2, int e, int w) { return ((((size_t)c * 2 + g2) * kNCls + e) * kNWS + w) * kNRepS; }   // g2: 0 C2, 1 C3
inline size_t SizeSWR() { return (size_t)kNSCfg * 2 * kNCls * kNWS * kNRepS; }
inline size_t IdxHq(int e, int p, int b) { return ((size_t)e * kNPP + p) * kNHq + b; }
inline size_t SizeHq() { return (size_t)kNCls * kNPP * kNHq; }

// ---------------------------------------------------------------------------------------------------------------
// 2. binary I/O
// ---------------------------------------------------------------------------------------------------------------
inline bool WriteArr(const std::string& f, const std::vector<double>& v) {
    FILE* fp = fopen(f.c_str(), "wb");
    if (!fp) { printf("bedge: cannot write %s\n", f.c_str()); return false; }
    const uint64_t n = v.size();
    fwrite(&n, sizeof(n), 1, fp);
    fwrite(v.data(), sizeof(double), v.size(), fp);
    fclose(fp);
    return true;
}
inline bool ReadArr(const std::string& f, std::vector<double>& v, size_t expect) {
    FILE* fp = fopen(f.c_str(), "rb");
    if (!fp) { printf("bedge: cannot read %s\n", f.c_str()); return false; }
    uint64_t n = 0;
    if (fread(&n, sizeof(n), 1, fp) != 1 || n != expect) { printf("bedge: %s has %llu values, expected %zu\n", f.c_str(), (unsigned long long)n, expect); fclose(fp); return false; }
    v.resize(n);
    const size_t got = fread(v.data(), sizeof(double), n, fp);
    fclose(fp);
    return got == n;
}

// ---------------------------------------------------------------------------------------------------------------
// 3. spectra
// ---------------------------------------------------------------------------------------------------------------
// CB rule options
struct CBOpt {
    double k = 1.0;   // CB x k
    int reb = 4;      // fine bins per CB cell (4 = 10 MeV; 2 = 5 MeV; 8 = 20 MeV)
    int mode = 0;     // 0: 2 sqrt(N++ N--); 1: N++ + N-- (like-sign sum)
};
// a scan::Spectra from counts us/pp/mm on nbin fine bins from lo (bin width bw), CB as ScanEngine::FillMany
inline scan::Spectra MakeSpec(const double* us, const double* pp, const double* mm, int nbin, double lo, double bw, const CBOpt& o = CBOpt()) {
    scan::Spectra sp;
    sp.nbin = nbin; sp.reb = o.reb; sp.nw = nbin / o.reb; sp.lo = lo; sp.hi = lo + nbin * bw; sp.bw = bw;
    const double* src[3] = {us, pp, mm};
    for (int k = 0; k < 3; ++k) {
        sp.n[k].assign(src[k], src[k] + nbin);
        sp.w[k].assign(sp.nw, 0.);
        for (int i = 0; i < nbin; ++i) if (i / o.reb < sp.nw) sp.w[k][i / o.reb] += sp.n[k][i];
        sp.nsel[k] = 0; for (int i = 0; i < nbin; ++i) sp.nsel[k] += sp.n[k][i];
    }
    sp.cbW.assign(sp.nw, 0.); sp.cbWVar.assign(sp.nw, 0.);
    for (int j = 0; j < sp.nw; ++j) {
        if (o.mode == 0) sp.cbW[j] = o.k * 2.0 * std::sqrt(sp.w[1][j] * sp.w[2][j]);
        else sp.cbW[j] = o.k * (sp.w[1][j] + sp.w[2][j]);
        sp.cbWVar[j] = o.k * o.k * (sp.w[1][j] + sp.w[2][j]);
    }
    sp.cb.assign(nbin, 0.); sp.cbVar.assign(nbin, 0.); sp.sig.assign(nbin, 0.); sp.sigVar.assign(nbin, 0.);
    for (int i = 0; i < nbin; ++i) {
        int j0, j1; double f; sp.Interp(i, j0, j1, f);
        const double a0 = (1 - f) / sp.reb, a1 = f / sp.reb;
        sp.cb[i] = a0 * sp.cbW[j0] + a1 * sp.cbW[j1];
        sp.cbVar[i] = a0 * a0 * sp.cbWVar[j0] + a1 * a1 * sp.cbWVar[j1];
        sp.sig[i] = sp.n[0][i] - sp.cb[i];
        sp.sigVar[i] = sp.n[0][i] + sp.cbVar[i];
    }
    return sp;
}
// a simulated spectrum (no CB): copy of wp6_fill.C SimSpectra
inline scan::Spectra SimSpec(const double* y, int nbin, double lo, double bw, int reb = 4) {
    scan::Spectra s;
    s.nbin = nbin; s.reb = reb; s.nw = nbin / reb; s.lo = lo; s.hi = lo + nbin * bw; s.bw = bw;
    for (int k = 0; k < 3; ++k) { s.n[k].assign(nbin, 0.); s.w[k].assign(s.nw, 0.); }
    for (int i = 0; i < nbin; ++i) s.n[0][i] = y[i];
    for (int i = 0; i < nbin; ++i) if (i / reb < s.nw) s.w[0][i / reb] += y[i];
    s.cbW.assign(s.nw, 0.); s.cbWVar.assign(s.nw, 0.);
    s.cb.assign(nbin, 0.); s.cbVar.assign(nbin, 0.);
    s.sig.assign(y, y + nbin); s.sigVar.assign(y, y + nbin);
    for (int k = 0; k < 3; ++k) s.nsel[k] = 0;
    for (int i = 0; i < nbin; ++i) s.nsel[0] += y[i];
    return s;
}
// spread-rule CB of a class spectrum: copy of wp6_fill.C SpreadCB (CB of the whole slice x the class share of the
// like-sign rows N++ + N-- in the 50 MeV band, bands from 0.30)
inline int BandOfX(double x) { return (int)std::floor((x - 0.30) / 0.05 + 1e-9); }
inline scan::Spectra SpreadCB(const scan::Spectra& all, const scan::Spectra& sl) {
    scan::Spectra s = sl;
    std::map<int, double> la, ls;
    for (int i = 0; i < sl.nbin; ++i) {
        const int b = BandOfX(sl.x(i));
        la[b] += all.n[1][i] + all.n[2][i];
        ls[b] += sl.n[1][i] + sl.n[2][i];
    }
    auto f = [&](double x) { const int b = BandOfX(x); return la[b] > 0 ? ls[b] / la[b] : 0.; };
    for (int i = 0; i < s.nbin; ++i) {
        const double fi = f(s.x(i));
        s.cb[i] = all.cb[i] * fi;
        s.cbVar[i] = all.cbVar[i] * fi * fi;
        s.sig[i] = s.n[0][i] - s.cb[i];
        s.sigVar[i] = s.n[0][i] + s.cbVar[i];
    }
    for (int j = 0; j < s.nw; ++j) {
        const double xc = s.lo + (j + 0.5) * s.reb * s.bw;
        const double fj = f(xc);
        s.cbW[j] = all.cbW[j] * fj;
        s.cbWVar[j] = all.cbWVar[j] * fj * fj;
    }
    return s;
}

// ---------------------------------------------------------------------------------------------------------------
// 4. edge-step fit
// ---------------------------------------------------------------------------------------------------------------
// step templates: 0 = scan::StepTemplate (erfc averaged uniformly over the slice, edge EdgeMeasured, q = 0);
// 1 = averaged also over the simulated M(ee) of the bump classes near the edge (edge E(m, q) of defs.h on the
//     scan scale); 2 = smeared "box in M4^2" (the phase-space shape at fixed M(pipi) for a massless gamma*:
//     density proportional to M4 between M4min and M4max, normalised to 1 at the edge), averaged over the slice
struct Tmpl {
    int type = 0;
    std::vector<double> q, wq;   // type 1: M(ee) values [GeV/c2] and weights
    std::map<long long, double> cache;   // value per bin centre (one Tmpl object per slice, smearing and scale)
};
inline double EdgeQ(double mp, double q, double scale, double scalePiPi) {
    const double me = scan::kMeta, mp0 = scan::kMpi0;
    const double mpr = mp / scalePiPi;
    const double E0 = (me * me + mp0 * mp0 - mpr * mpr) / (2. * me);
    const double p2 = E0 * E0 - mp0 * mp0;
    if (!(p2 >= 0.)) return -1.;
    const double a = me * me - me * (E0 - std::sqrt(p2)) * (1. - q * q / (mp0 * mp0));
    return a > 0. ? scale * std::sqrt(a) : -1.;
}
inline double EdgeLo(double mp, double scale, double scalePiPi) {
    const double me = scan::kMeta, mp0 = scan::kMpi0;
    const double mpr = mp / scalePiPi;
    const double E0 = (me * me + mp0 * mp0 - mpr * mpr) / (2. * me);
    const double p2 = E0 * E0 - mp0 * mp0;
    if (!(p2 >= 0.)) return -1.;
    const double a = me * me - me * (E0 + std::sqrt(p2));
    return a > 0. ? scale * std::sqrt(a) : -1.;
}
inline double TmplValue(double m, const scan::StepCfg& c, Tmpl& t) {
    if (t.type == 0 || c.useEdge) return scan::StepTemplate(m, c);
    const long long key = std::llround(m * 1e8);
    auto it = t.cache.find(key);
    if (it != t.cache.end()) return it->second;
    const double qq = 1.0 / (std::sqrt(2.0) * c.s);
    const int n = 21;
    double s = 0, wsum = 0;
    for (int i = 0; i < n; ++i) {
        const double mp = c.mpLo + (i + 0.5) * (c.mpHi - c.mpLo) / n;
        if (t.type == 1) {
            for (size_t j = 0; j < t.q.size(); ++j) {
                const double e = EdgeQ(mp, t.q[j], c.scale, c.scalePiPi);
                if (e <= 0) continue;
                s += t.wq[j] * 0.5 * std::erfc((m - e) * qq); wsum += t.wq[j];
            }
        } else {
            // smeared box: int_{Mlo}^{Mhi} (M/Mhi) G(m - M; s) dM, by 200-point midpoint rule
            const double Mhi = scan::EdgeMeasured(mp, c.scale, c.scalePiPi), Mlo = EdgeLo(mp, c.scale, c.scalePiPi);
            if (Mhi <= 0 || Mlo <= 0) continue;
            const int nq = 200; double v = 0;
            const double a = std::max(Mlo, m - 6 * c.s), b = std::min(Mhi, m + 6 * c.s);
            if (b > a) {
                const double d = (b - a) / nq;
                for (int k = 0; k < nq; ++k) {
                    const double M = a + (k + 0.5) * d;
                    v += (M / Mhi) * std::exp(-0.5 * (m - M) * (m - M) / (c.s * c.s)) * d;
                }
                v /= std::sqrt(2 * M_PI) * c.s;
            }
            s += v; wsum += 1;
        }
    }
    const double val = wsum > 0 ? s / wsum : 0.;
    t.cache[key] = val;
    return val;
}
// the fit: copy of scan::FitEdgeStep; Vfix (npts x npts, optional) replaces the iterated covariance; Vout gets
// the final covariance; tpl selects the step template
inline scan::StepFit FitStep(const scan::Spectra& sp, const scan::StepCfg& c, Tmpl* tpl = nullptr,
                             const std::vector<double>* Vfix = nullptr, std::vector<double>* Vout = nullptr) {
    using namespace scan;
    StepFit r; r.cfg = c;
    for (int k = 0; k < kMaxPar; ++k) r.par[k] = 0;
    if (c.useEdge) { r.edge = r.edgeLo = r.edgeHi = c.edge; }
    else {
        r.edge = EdgeMeasured(0.5 * (c.mpLo + c.mpHi), c.scale, c.scalePiPi);
        r.edgeHi = EdgeMeasured(c.mpLo, c.scale, c.scalePiPi);
        r.edgeLo = EdgeMeasured(c.mpHi, c.scale, c.scalePiPi);
        if (r.edge <= 0) return r;
    }
    r.lo = c.lo > 0 ? c.lo : r.edge - c.halfRange;
    r.hi = c.hi > 0 ? c.hi : r.edge + c.halfRange;
    if (c.withEta >= 0) r.withEta = (c.withEta == 1);
    else if (r.hi <= c.etaMean - 3.0 * c.etaSigma) r.withEta = false;
    else if (c.hi > 0) r.withEta = true;
    else if (r.hi < c.etaMean - c.etaSigma) { r.hi = c.etaMean - 3.0 * c.etaSigma; r.withEta = false; }
    else { r.hi = std::max(r.hi, c.etaMean + 3.0 * c.etaSigma); r.withEta = true; }
    r.lo = sp.lo + sp.bw * std::floor((r.lo - sp.lo) / sp.bw + 0.5);
    r.hi = sp.lo + sp.bw * std::floor((r.hi - sp.lo) / sp.bw + 0.5);
    r.x0 = 0.5 * (r.lo + r.hi); r.dx = 0.5 * (r.hi - r.lo);
    std::vector<int> idx;
    for (int i = 0; i < sp.nbin; ++i) if (sp.x(i) >= r.lo && sp.x(i) < r.hi) idx.push_back(i);
    const int npts = (int)idx.size();
    const int npol = c.order + 1;
    const int np = npol + (r.withEta ? 1 : 0) + 1;
    r.np = np;
    if (np > kMaxPar || npts <= np) return r;
    if (Vfix && (int)Vfix->size() != npts * npts) return r;
    std::vector<double> F(npts * np), F0(npts * (np - 1)), y(npts), var(npts), ad(npts), V0(npts * npts), V;
    for (int q = 0; q < npts; ++q) {
        const int i = idx[q];
        const double xc = sp.x(i), t = (xc - r.x0) / r.dx;
        double tp = 1; int k = 0;
        for (int j = 0; j < npol; ++j) { F[q * np + k++] = tp; tp *= t; }
        if (r.withEta) F[q * np + k++] = GaussFrac(xc - 0.5 * sp.bw, xc + 0.5 * sp.bw, c.etaMean, c.etaSigma);
        const double st = tpl ? TmplValue(xc, c, *tpl) : StepTemplate(xc, c);
        F[q * np + k++] = st;
        for (int j = 0; j < np - 1; ++j) F0[q * (np - 1) + j] = F[q * np + j];
        r.tx.push_back(xc); r.ty.push_back(st);
        y[q] = sp.sig[i]; ad[q] = sp.cb[i];
        var[q] = std::max(y[q] + ad[q], 1.0);
    }
    LinRes a, b;
    if (Vfix) {
        V = *Vfix; a = LinFitCov(npts, np, y.data(), V, F.data());
        if (!a.ok) return r;
        V = *Vfix; b = LinFitCov(npts, np - 1, y.data(), V, F0.data());
        if (!b.ok) return r;
        if (Vout) *Vout = *Vfix;
    } else {
        for (int q = 0; q < npts; ++q)
            for (int p = 0; p <= q; ++p) V0[q * npts + p] = V0[p * npts + q] = sp.CbCov(idx[q], idx[p]);
        auto solve = [&](int npar, const double* FF) -> LinRes {
            V = V0;
            for (int q = 0; q < npts; ++q) V[q * npts + q] += var[q];
            return LinFitCov(npts, npar, y.data(), V, FF);
        };
        a = solve(np, F.data());
        if (!a.ok) return r;
        for (int it = 0; it < 6; ++it) {
            const double hOld = a.par[np - 1];
            for (int q = 0; q < npts; ++q) {
                double f = 0;
                for (int k = 0; k < np; ++k) f += a.par[k] * F[q * np + k];
                var[q] = std::max(f + ad[q], 1.0);
            }
            a = solve(np, F.data());
            if (!a.ok) return r;
            if (std::fabs(a.par[np - 1] - hOld) < 1e-4 * std::sqrt(a.cov[np - 1][np - 1])) break;
        }
        b = solve(np - 1, F0.data());
        if (!b.ok) return r;
        if (Vout) { *Vout = V0; for (int q = 0; q < npts; ++q) (*Vout)[q * npts + q] += var[q]; }
    }
    for (int k = 0; k < np; ++k) r.par[k] = a.par[k];
    for (int q = 0; q < npts; ++q) {
        double f = 0;
        for (int k = 0; k < np; ++k) f += a.par[k] * F[q * np + k];
        r.fy.push_back(f);
        r.fy0.push_back(f - a.par[np - 1] * F[q * np + np - 1]);
    }
    r.h = a.par[np - 1]; r.eh = std::sqrt(a.cov[np - 1][np - 1]);
    r.chi2 = a.chi2; r.chi2NoStep = b.chi2; r.ndf = npts - np;
    r.dchi2 = b.chi2 - a.chi2;
    r.bkgAtEdge = r.Poly(r.edge);
    if (r.withEta) { r.Seta = a.par[npol]; r.eSeta = std::sqrt(a.cov[npol][npol]); }
    r.ok = true;
    return r;
}

}  // namespace bedge

#endif
