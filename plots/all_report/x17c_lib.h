// ============================================================================
// x17c_lib.h  --  WP9 X-c (prefix x17c): the FROZEN X17 procedure as code (inputs, configuration, side-band fit,
//                 masked-row prediction, limits). Used by x17c_freeze.C, x17c_toys.C, x17c_summary.C and, after the
//                 freeze, by X-d (observed limit). Written 2026-11-14. Model: x17c_model.h (= x17_model.h of x17ab plus
//                 the normalisation switch, line-shape nuisances and the variance factor).
// ============================================================================
// Blinding (D19): the data cells of X2D with M(ee) 12-28 MeV/c2 (rows 6-13, 1-based; the 2 MeV bins that overlap 13-27
// MeV) are never read here: LoadSideBands() skips them and only checks that hists_v1.root holds the mask value -1 there.
// ============================================================================
#ifndef ALLREPORT_X17C_LIB_H
#define ALLREPORT_X17C_LIB_H

#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH2.h"
#include "TH3.h"

#include "defs.h"
#include "hists_io.h"
#include "x17c_model.h"

namespace x17cl {

const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/x17c";
const char* const kFrozenDir = "/hdd3/przygoda/HADES/FAT_out/all_report/x17";
const char* const kFig = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
const char* const kTab = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all";
const char* const kWp5Results = "/hdd3/przygoda/HADES/FAT_out/all_report/wp5/results.txt";
const char* const kConfigFile = "/hdd3/przygoda/HADES/FAT_out/all_report/x17c/x17c_config.txt";
const double kBR = x17c::kBRdir;   // 2.68e-4
const int kMaskLo = 5, kMaskHi = 12;   // 0-based X2D rows 12-28 MeV/c2 (never read)

// mass grid of the scan (X-d): 10-40 MeV/c2 in 1 MeV/c2, 0.5 MeV/c2 in 15-20, plus the quoted masses 16.7 and 16.85
inline std::vector<double> MassGrid() {
    std::vector<double> m;
    for (int k = 10; k <= 14; ++k) m.push_back(k);
    for (int k = 0; k <= 10; ++k) m.push_back(15. + 0.5 * k);
    for (int k = 21; k <= 40; ++k) m.push_back(k);
    m.push_back(16.7); m.push_back(16.85);
    std::sort(m.begin(), m.end());
    return m;
}
const double kQuoted[3] = {16.7, 16.85, 17.0};

// ---------------------------------------------------------------- WP5 inputs (results.txt)
struct Wp5 {
    bool ok = false;
    double a = 0, aStat = 0, aTot = 0, aBias = 0, g = 0, gTot = 0, kB = 0, kBTot = 0, b = 0, bTot = 0;
    double chi2scale = 1.61;   // sqrt(GLS chi2/ndf) of F-main (wp5_results.md)
    double Ypeak = 0, YpeakTot = 0, YpeakQuad = 0;
    std::vector<double> B[2];   // background per X2D M4 bin, E1 (prediction) and E2
};
inline Wp5 ReadWp5() {
    Wp5 w;
    FILE* f = fopen(kWp5Results, "r");
    if (!f) { printf("x17cl: cannot open %s\n", kWp5Results); return w; }
    char line[4096];
    w.B[0].assign(alldefs::kX2Dnx, 0.); w.B[1].assign(alldefs::kX2Dnx, 0.);
    int nb = 0;
    while (fgets(line, sizeof line, f)) {
        double a, b, c, d;
        int e, i;
        if (sscanf(line, "x17_aeta %lf stat %lf total %lf bias_corr %lf", &a, &b, &c, &d) == 4) { w.a = a; w.aStat = b; w.aTot = c; w.aBias = d; }
        else if (sscanf(line, "x17_g %lf total %lf", &a, &b) == 2) { w.g = a; w.gTot = b; }
        else if (sscanf(line, "x17_kappaB %lf total %lf", &a, &b) == 2) { w.kB = a; w.kBTot = b; }
        else if (sscanf(line, "x17_b %lf total %lf", &a, &b) == 2) { w.b = a; w.bTot = b; }
        else if (sscanf(line, "x17bkg E%d %d %lf %lf %lf", &e, &i, &a, &b, &c) == 5) {
            if (e >= 1 && e <= 2 && i >= 1 && i <= alldefs::kX2Dnx) { w.B[e - 1][i - 1] = c; ++nb; }
        } else if (strncmp(line, "derived FMAIN Y_peak ", 21) == 0) {
            double v, tot;
            const char* p = strstr(line, "total ");
            if (sscanf(line + 21, "%lf", &v) == 1 && p && sscanf(p + 6, "%lf", &tot) == 1) { w.Ypeak = v; w.YpeakTot = tot; }
        } else if (strncmp(line, "derived FMAIN_quad Y_peak ", 26) == 0) {
            double v;
            if (sscanf(line + 26, "%lf", &v) == 1) w.YpeakQuad = v;
        }
    }
    fclose(f);
    w.ok = w.a > 0 && w.g > 0 && nb == 2 * alldefs::kX2Dnx;
    if (!w.ok) printf("x17cl: WP5 inputs incomplete (a %.3f g %.3f, %d background bins)\n", w.a, w.g, nb);
    return w;
}

// ---------------------------------------------------------------- data (side bands only)
// n: [nx*ny] side-band rows filled, masked rows 0; use: 1 for side-band cells. Returns false if the mask is not intact.
inline bool LoadSideBands(TFile* f, std::vector<double>& n, std::vector<char>& use, const char* sel = "S0") {
    TH2* h = (TH2*)f->Get(Form("data/%s/all/N_X2D", sel));
    if (!h) { printf("x17cl: data N_X2D missing\n"); return false; }
    const int nx = h->GetNbinsX(), ny = h->GetNbinsY();
    n.assign((size_t)nx * ny, 0.);
    use.assign((size_t)nx * ny, 0);
    int bad = 0;
    for (int j = 0; j < ny; ++j) {
        const bool masked = (j >= kMaskLo && j <= kMaskHi);
        for (int i = 0; i < nx; ++i) {
            if (masked) {   // only the mask value is checked; a real count is never copied or printed
                if (h->GetBinContent(i + 1, j + 1) != -1.) ++bad;
                continue;
            }
            n[(size_t)j * nx + i] = h->GetBinContent(i + 1, j + 1);
            use[(size_t)j * nx + i] = 1;
        }
    }
    if (bad) { printf("x17cl: BLINDING MASK NOT INTACT in %d cells -- stop\n", bad); return false; }
    return true;
}

// ---------------------------------------------------------------- line of mass m with the line-shape derivatives
const double kWidthSig = 0.10, kScaleSig = 0.02;   // spec WP9: resolution +-10 %, M(ee) scale +-2 %
inline bool SetLine(x17c::Inputs& I, double m, int rxVar = 0, double widthFac = 1., double scaleFac = 1.) {
    double rx, d;
    std::vector<double> Xp, Xm, Sp, Sm;
    if (!x17c::LineFor(m, I.G, I.X, rx, rxVar, scaleFac, widthFac - 1.)) return false;
    I.RX = rx;
    x17c::LineFor(m, I.G, Xp, d, rxVar, scaleFac, widthFac * (1 + kWidthSig) - 1.);
    x17c::LineFor(m, I.G, Xm, d, rxVar, scaleFac, widthFac * (1 - kWidthSig) - 1.);
    x17c::LineFor(m, I.G, Sp, d, rxVar, scaleFac * (1 + kScaleSig), widthFac - 1.);
    x17c::LineFor(m, I.G, Sm, d, rxVar, scaleFac * (1 - kScaleSig), widthFac - 1.);
    I.dXw.assign(I.G.ny, 0.); I.dXs.assign(I.G.ny, 0.);
    for (int j = 0; j < I.G.ny; ++j) { I.dXw[j] = 0.5 * (Xp[j] - Xm[j]); I.dXs[j] = 0.5 * (Sp[j] - Sm[j]); }
    return true;
}

// ---------------------------------------------------------------- inputs
struct BuildOpt {
    int bshape = 0;          // 0 = WP5 background shapes (NOMINAL); 1 = simulated G_D + G_C exp(pol2) (x17ab); 2 = exp(pol1);
                             // 3 = exp(pol3); 4 = simulated with kappa_B (WP5) ; (bkgMode 2 = flat is a Config choice)
    int nbcl = 2;
    int cbMode = 0;          // 0 spread rule, 1 factorised, 2 per-row geometric mean
    int cbRep = -1;          // CB replica
    int simRep = -1;         // C4 replica
    std::string var = "nom"; // simulation variant directory for the peak/direct templates
};
// rxHist (optional): R_X of the histograms (C6 / C5 of the variant directory), to carry a variant's R_X change onto the line
inline bool Build(TFile* f, x17c::Inputs& I, const Wp5& w, const BuildOpt& b, double aStart, double* rxHist = nullptr) {
    x17c::HistOpt o;
    o.aEta = 1.;
    o.nbcl = b.nbcl;
    o.cbMode = b.cbMode;
    o.cbRep = b.cbRep;
    o.simRep = b.simRep;
    o.smoothOrder = (b.bshape == 2) ? 1 : (b.bshape == 3 ? 3 : 2);
    if (b.bshape == 4) o.kappaB = w.kB;
    if (!x17c::BuildFromHists(f, I, o, alldefs::kNbToRaw)) return false;
    if (b.var != "nom") {   // peak/direct templates, D and R_X of a simulation variant; everything else nominal
        x17c::HistOpt ov = o; ov.var = b.var;
        x17c::Inputs Iv;
        if (!x17c::BuildFromHists(f, Iv, ov, alldefs::kNbToRaw)) return false;
        I.D = Iv.D; I.gP = Iv.gP; I.gX = Iv.gX; I.RX = Iv.RX; I.N5start = Iv.N5start;
    }
    I.N5ref = I.N5start;   // direct rows (C5, V3), full M(ee), X2D M4 range, at a_eta = 1
    if (rxHist) *rxHist = I.RX;
    if (b.bshape == 0 && b.nbcl == 2) {
        for (int c = 0; c < 2; ++c) { I.B[c] = w.B[c]; x17c::Normalise(I.B[c]); }
    }
    I.N5start = aStart * I.N5ref;
    I.N4start *= aStart;
    return true;
}
// simulated G_D + G_C per M(ee) row (raw rows): reference for the masked-row prediction
inline std::vector<double> SimRow(TFile* f, const x17c::Inputs& I) {
    std::vector<double> s(I.G.ny, 0.);
    TH2* gd = (TH2*)f->Get("sim/S0__nom/GD_X2D");
    TH2* gc = (TH2*)f->Get("sim/S0__nom/Gconv_X2D");
    for (int j = 0; j < I.G.ny; ++j) for (int i = 1; i <= I.G.nx; ++i) s[j] += (gd->GetBinContent(i, j + 1) + gc->GetBinContent(i, j + 1)) * alldefs::kNbToRaw;
    return s;
}
// relative simulation error of the G_D + G_C row sums (weights of the masked-row prediction)
inline std::vector<double> SimRowRelErr(TFile* f, const x17c::Inputs& I) {
    std::vector<double> s(I.G.ny, 0.), e2(I.G.ny, 0.);
    TH2* gd = (TH2*)f->Get("sim/S0__nom/GD_X2D");
    TH2* gc = (TH2*)f->Get("sim/S0__nom/Gconv_X2D");
    for (int j = 0; j < I.G.ny; ++j) for (int i = 1; i <= I.G.nx; ++i) {
        s[j] += gd->GetBinContent(i, j + 1) + gc->GetBinContent(i, j + 1);
        e2[j] += std::pow(gd->GetBinError(i, j + 1), 2) + std::pow(gc->GetBinError(i, j + 1), 2);
    }
    for (int j = 0; j < I.G.ny; ++j) s[j] = s[j] > 0 ? std::sqrt(e2[j]) / s[j] : 0.;
    return s;
}
// weighted least-squares fit of ln(s_j) by a polynomial of order 'order' in y = (M(ee) - 30 MeV)/30 MeV over the rows
// jlo..jhi; returns the curve on all rows (chi2/ndf of the fit in *chi2ndf)
inline std::vector<double> SmoothLnPol(const std::vector<double>& s, const std::vector<double>& e2, const x17c::Inputs& I, int order, int jlo, int jhi,
                                       double* chi2ndf) {
    const int nc = order + 1;
    std::vector<double> A((size_t)nc * nc, 0.), b(nc, 0.);
    int nfit = 0;
    for (int j = jlo; j <= jhi && j < I.G.ny; ++j) {
        if (!(s[j] > 0) || !(e2[j] > 0)) continue;
        const double y = (I.G.yc(j) - 0.030) / 0.030, z = std::log(s[j]), w = s[j] * s[j] / e2[j];
        std::vector<double> xp(nc, 1.);
        for (int q = 1; q < nc; ++q) xp[q] = xp[q - 1] * y;
        for (int u = 0; u < nc; ++u) { b[u] += w * xp[u] * z; for (int v = 0; v < nc; ++v) A[(size_t)u * nc + v] += w * xp[u] * xp[v]; }
        ++nfit;
    }
    if (nfit <= nc || !x17c::Model::Chol(A, b, nc)) return s;
    std::vector<double> out(I.G.ny, 0.);
    double c2 = 0;
    for (int j = 0; j < I.G.ny; ++j) {
        const double y = (I.G.yc(j) - 0.030) / 0.030;
        double z = 0, xp = 1;
        for (int q = 0; q < nc; ++q) { z += b[q] * xp; xp *= y; }
        out[j] = std::exp(z);
        if (j >= jlo && j <= jhi && e2[j] > 0 && s[j] > 0) c2 += std::pow(s[j] - out[j], 2) / e2[j];
    }
    if (chi2ndf) *chi2ndf = c2 / std::max(1, nfit - nc);
    return out;
}
// LOCALLY SMOOTHED simulated G_D + G_C per M(ee) row (variant reference of the masked-row prediction): for M(ee) >= 8 MeV/c2
// G_D = exp(pol2) fitted over 8-62 MeV/c2 and G_C = exp(pol1) fitted over 8-28 MeV/c2 (G_C is ~0 above); below 8 MeV/c2
// the raw sums (3-5 % errors). chi2D/chi2C: chi2/ndf of the two fits.
inline std::vector<double> SimRowSmooth(TFile* f, const x17c::Inputs& I, double* chi2D = nullptr, double* chi2C = nullptr) {
    std::vector<double> sD(I.G.ny, 0.), eD(I.G.ny, 0.), sC(I.G.ny, 0.), eC(I.G.ny, 0.);
    TH2* gd = (TH2*)f->Get("sim/S0__nom/GD_X2D");
    TH2* gc = (TH2*)f->Get("sim/S0__nom/Gconv_X2D");
    const double k = alldefs::kNbToRaw;
    for (int j = 0; j < I.G.ny; ++j) for (int i = 1; i <= I.G.nx; ++i) {
        sD[j] += gd->GetBinContent(i, j + 1) * k; eD[j] += std::pow(gd->GetBinError(i, j + 1) * k, 2);
        sC[j] += gc->GetBinContent(i, j + 1) * k; eC[j] += std::pow(gc->GetBinError(i, j + 1) * k, 2);
    }
    int j8 = 0, j28 = 0;
    for (int j = 0; j < I.G.ny; ++j) { if (I.G.ye[j] < 0.008 - 1e-9) j8 = j + 1; if (I.G.ye[j + 1] <= 0.028 + 1e-9) j28 = j; }
    const std::vector<double> mD = SmoothLnPol(sD, eD, I, 2, j8, I.G.ny - 1, chi2D);
    const std::vector<double> mC = SmoothLnPol(sC, eC, I, 1, j8, j28, chi2C);
    std::vector<double> out(I.G.ny, 0.);
    for (int j = 0; j < I.G.ny; ++j) out[j] = j < j8 ? sD[j] + sC[j] : mD[j] + mC[j];
    return out;
}
// like-sign rows (N++ + N--) per M(ee) row (variant reference)
inline std::vector<double> LsRow(TFile* f, const x17c::Inputs& I) {
    std::vector<double> s(I.G.ny, 0.);
    TH2* pp = (TH2*)f->Get("data/S0/all/Npp_X2D");
    TH2* mm = (TH2*)f->Get("data/S0/all/Nmm_X2D");
    for (int j = 0; j < I.G.ny; ++j) for (int i = 1; i <= I.G.nx; ++i) s[j] += pp->GetBinContent(i, j + 1) + mm->GetBinContent(i, j + 1);
    return s;
}

// ---------------------------------------------------------------- configuration
struct Frozen {
    // constraints (WP5) and nuisance sizes; filled by x17c_freeze.C and stored in x17c_config.txt
    double a0 = 0.382, aSig = 0.164, g0 = 1.717, gSig = 0.645;
    double kE1 = 0.5, kE2 = 0.05;
    double cTail = 0.194;
    double xSig1 = 0.20, xSig2 = 0.45, xSig3 = 0.30;   // X17 yield nuisance (log-normal), normalisations (1) / (2) / (3)
    double qScale = 1.;                  // variance factor applied to q~
    double aP = 0, aPSig = 0;            // peak-based eta scale for normalisation (3)
    int bkgMode = 4;
    int bshape = 0;
    int fillOrder = 2;
};
// norm 1: X17 rows = r N5 R_X (N5 fitted); 2: r a_eta(WP5) N5ref R_X; 3: r a_P N5ref R_X. The a_eta parameter of the
// eta-related background terms has the WP5 prior in all three.
inline x17c::Config MakeConfig(const Frozen& z, int norm) {
    x17c::Config C;
    C.kSig[0] = z.kE1; C.kSig[1] = z.kE2;
    C.a0 = z.a0;
    C.aSig = z.aSig;
    C.gC = z.g0; C.gSig = z.gSig;
    C.bkgMode = z.bkgMode;
    C.cTailSig = z.cTail;
    C.xSig = norm == 1 ? z.xSig1 : (norm == 2 ? z.xSig2 : z.xSig3);
    C.wSig = 1.; C.sSig = 1.;
    C.qScale = z.qScale;
    C.xNorm = norm == 1 ? 0 : 1;
    C.xA = norm == 3 ? z.aP : z.a0;
    return C;
}
inline bool WriteConfig(const Frozen& z, const char* file = kConfigFile) {
    FILE* f = fopen(file, "w");
    if (!f) return false;
    fprintf(f, "a0 %.6g\naSig %.6g\ng0 %.6g\ngSig %.6g\nkE1 %.6g\nkE2 %.6g\ncTail %.6g\nxSig1 %.6g\nxSig2 %.6g\nxSig3 %.6g\nqScale %.6g\naP %.6g\naPSig %.6g\nbkgMode %d\nbshape %d\nfillOrder %d\n",
            z.a0, z.aSig, z.g0, z.gSig, z.kE1, z.kE2, z.cTail, z.xSig1, z.xSig2, z.xSig3, z.qScale, z.aP, z.aPSig, z.bkgMode, z.bshape, z.fillOrder);
    fclose(f);
    return true;
}
inline bool ReadConfig(Frozen& z, const char* file = kConfigFile) {
    FILE* f = fopen(file, "r");
    if (!f) { printf("x17cl: cannot read %s\n", file); return false; }
    char k[64];
    double v;
    std::map<std::string, double> m;
    while (fscanf(f, "%63s %lf", k, &v) == 2) m[k] = v;
    fclose(f);
    auto g = [&](const char* n, double& x) { if (m.count(n)) x = m[n]; };
    g("a0", z.a0); g("aSig", z.aSig); g("g0", z.g0); g("gSig", z.gSig); g("kE1", z.kE1); g("kE2", z.kE2); g("cTail", z.cTail);
    g("xSig1", z.xSig1); g("xSig2", z.xSig2); g("xSig3", z.xSig3); g("qScale", z.qScale); g("aP", z.aP); g("aPSig", z.aPSig);
    double t;
    if (m.count("bkgMode")) { t = m["bkgMode"]; z.bkgMode = (int)t; }
    if (m.count("bshape")) { t = m["bshape"]; z.bshape = (int)t; }
    if (m.count("fillOrder")) { t = m["fillOrder"]; z.fillOrder = (int)t; }
    return true;
}

// ---------------------------------------------------------------- side-band fit and the background-only expectation
struct SideBand {
    x17c::FitResult fit;          // r = 0, side-band cells only
    std::vector<double> pFill;    // parameters with the masked n_j predicted
    double chi2 = 0;              // Pearson chi2 of the side-band cells
    int ndf = 0, nfree = 0, ncell = 0;
    double dev = 0;               // deviance (-2 ln L without the penalty) of the side-band cells
    double chi2c[2] = {0, 0};     // per class
    int nc[2] = {0, 0};
    std::vector<double> chi2row;  // per M(ee) row
    int nFilled = 0;
};
// refRel (optional): relative uncertainty of the reference rows, added in quadrature to the relative error of n_j in the
// weights of the masked-row prediction
inline SideBand FitSideBands(const x17c::Inputs& I, const x17c::Config& Csb, const std::vector<double>& n, const std::vector<double>& ref, int fillOrder,
                             const std::vector<double>* refRel = nullptr) {
    SideBand S;
    x17c::Model M(I, Csb);
    const x17c::Glob g = x17c::GlobOf(Csb);
    S.fit = M.Fit(n, 0., false, g);
    if (!S.fit.ok) {   // second attempt from the result
        x17c::FitResult f2 = M.Fit(n, 0., false, g, &S.fit.p, 600);
        if (f2.F <= S.fit.F) S.fit = f2;
    }
    const std::vector<double>& p = S.fit.p;
    S.chi2row.assign(I.G.ny, 0.);
    for (int j = 0; j < I.G.ny; ++j)
        for (int i = 0; i < I.G.nx; ++i) {
            const int b = I.G.idx(i, j);
            if (!M.Used(b)) continue;
            const double mu = M.Mu(p, i, j), nn = n[b];
            if (!(mu > 0)) continue;
            const double c = (nn - mu) * (nn - mu) / mu;
            S.chi2 += c; S.chi2row[j] += c; S.chi2c[I.kcl[j]] += c; S.nc[I.kcl[j]]++;
            S.dev += 2 * (mu - nn + (nn > 0 ? nn * std::log(nn / mu) : 0.));
            ++S.ncell;
        }
    for (size_t k = 0; k < S.fit.err.size(); ++k) if (S.fit.err[k] > 0) ++S.nfree;
    S.ndf = S.ncell - S.nfree;
    S.pFill = p;
    std::vector<double> nErr(I.G.ny, 0.);
    for (int j = 0; j < I.G.ny; ++j) {
        nErr[j] = S.fit.err[j];
        if (refRel && p[j] > 0 && nErr[j] > 0) nErr[j] = p[j] * std::hypot(nErr[j] / p[j], (*refRel)[j]);
    }
    S.nFilled = x17c::FillMaskedRows(M, S.pFill, ref, fillOrder, &nErr);
    S.pFill[M.L.iR] = 0.;
    return S;
}

// ---------------------------------------------------------------- content of the counting window X17c (spec Sect. 2)
// M4 0.5312-0.5566 GeV/c2 (data scale) x M(ee) 15.5-24.5 MeV/c2; X2D cells taken with their overlap fraction (uniform
// density inside a 2.5 MeV x 2 MeV cell assumed)
const double kX17cM4lo = 0.5312, kX17cM4hi = 0.5566, kX17cMeelo = 0.0155, kX17cMeehi = 0.0245;
inline double WindowX17c(const x17c::Grid& G, const std::vector<double>& v) {
    double s = 0;
    for (int j = 0; j < G.ny; ++j) {
        const double fy = std::max(0., std::min(G.ye[j + 1], kX17cMeehi) - std::max(G.ye[j], kX17cMeelo)) / (G.ye[j + 1] - G.ye[j]);
        if (fy <= 0) continue;
        for (int i = 0; i < G.nx; ++i) {
            const double fx = std::max(0., std::min(G.xe[i + 1], kX17cM4hi) - std::max(G.xe[i], kX17cM4lo)) / (G.xe[i + 1] - G.xe[i]);
            if (fx > 0) s += fx * fy * v[G.idx(i, j)];
        }
    }
    return s;
}

// ---------------------------------------------------------------- limit of one data set (all cells), in BR units
struct Lim {
    bool ok = false;
    double ul = 0, exp[5] = {0, 0, 0, 0, 0}, rhat = 0, rerr = 0, N5 = 0, nX = 0, fac = 0;   // fac: BR per unit r
};
// g: centres of the Gaussian constraints (global observables); nullptr = the nominal centres of C (data, Asimov)
inline Lim LimitBR(const x17c::Inputs& I, const x17c::Config& C, const std::vector<double>& n, bool doExp, double scale = 0.1,
                   const x17c::Glob* g = nullptr) {
    x17c::Model M(I, C);
    x17c::Limit L = x17c::LimitAsym(M, n, g ? *g : x17c::GlobOf(C), scale, doExp);
    Lim r;
    r.ok = L.ok;
    r.fac = kBR;
    r.ul = L.ul * kBR;
    for (int k = 0; k < 5; ++k) r.exp[k] = L.exp[k] * kBR;
    r.rhat = L.rhat * kBR; r.rerr = L.rhatErr * kBR;
    if (!L.prof.hat.p.empty()) {
        r.N5 = L.prof.hat.p[M.L.iN5];
        r.nX = L.ul * M.Nnorm(L.prof.zero.p) * I.RX;   // X17 rows at the limit (at the r = 0 fit's yield factor)
    }
    return r;
}

}  // namespace x17cl

#endif  // ALLREPORT_X17C_LIB_H
