// ============================================================================
// opt_sculpt.C  --  sculpting checks of the cuts on the data: can a cut make a
//                   peak (or a step) out of a smooth background?
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
//
// What it does (method of plots/scan_peak.C, "guard against sculpting", here on
// the base selection pre = cor_ok && isBest >= 0 && vz > -160 mm)
//   1. three row sets: US = unlike-sign rows, LS = like-sign rows (e+e+ and e-e-
//      added), MIX = mixed rows: the pion pair of one unlike-sign row with the
//      lepton pair of a row of ANOTHER event group (20 fixed offsets in the row
//      order; M4, oa_lab, the eta-rest angle and MM4 recomputed from the corrected
//      momenta; M(pipi), dz_pi, the pion vertex flags and vz from the pion row,
//      OA and M(ee) from the lepton row). Mixed rows have no eta and no pair
//      correlation, but the kinematic limits hold (M4 > M(pipi) + M(ee),
//      MM4 < sqrt(s) - M4).
//   2. pass fraction f(M4) of every scanned cut value and of the named selections
//      for the three sets; turn-on interval (M10-M90) and the change of f inside
//      the fit range 0.5225-0.6100 (dfit); flag G = like-sign turn-on inside the
//      fit range with dfit >= 0.05;
//   3. fake peak: y(M4) = Bref(M4) x f_MIX(M4), Bref = smooth background of pre
//      (pol2 of the Poisson fit of the unsubtracted N), fitted with the frozen
//      Gaussian + pol2. A pol2 background stays pol2 only if f is flat or linear in
//      the fit range; S_fake = Gaussian content [rows at the data normalisation];
//      error = larger of the binomial error and half the difference between the
//      two halves of the mixing offsets; flag M = |S_fake| > 0.25 sigma_S(pre);
//   4. like-sign null test of the named selections (Poisson fit of the N++ and
//      N-- spectra with the frozen Gaussian + pol2: expected 0);
//   5. figures: pass fractions against M4 for four cuts (M(pipi) < 0.42, MM4 <
//      2.85, eta-rest angle > 140 deg, oa_lab < 55) and S_fake against the cut value.
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root (scan::kCompact)
// Outputs: /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/opt_sculpt.txt
//          /home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/opt_sculpt_passfrac.pdf, opt_sculpt_fake.pdf (+ .png)
//
// Run (from the repository root /home/przygoda/HADES/FAT/FAT), about 2 min:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_opt", kTRUE)' 'plots/exp_report/opt_sculpt.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "opt_common.h"

namespace optsc {

using namespace opt;

// ---- cuts tested ---------------------------------------------------------------
struct Test { std::string fam, lab; Cuts c; double x = 0; };

// ---- mixed rows (copy of plots/scan_peak.C, pool restricted to the base) ----------
const int kNOff = 20;
const int kOff[kNOff] = {1009, 2003, 3001, 4001, 5003, 6007, 7001, 8009, 9001, 10007,
                         11003, 12007, 13001, 14009, 15013, 16001, 17011, 18013, 19001, 20011};
struct MixPool {
    const std::vector<Row>* rows = nullptr;
    std::vector<int> idx;
    std::vector<float> pp, ee;
    double pb = 0, eIni = 0;
};
struct MixKin { double m4 = 0, oal = 0, mm = 0, oae = 0; };
void BuildPool(const scan::Data& d, const Cuts& base, MixPool& mp) {
    const double mPi = 1e3 * scan::kMpi, mE = 1e3 * scan::kMe, mP = 1e3 * scan::kMp, d2r = TMath::DegToRad();
    const double eb = scan::kTbeam + mP;
    mp.pb = std::sqrt(eb * eb - mP * mP); mp.eIni = eb + mP;
    mp.rows = &d.s[0];
    const std::vector<Row>& v = d.s[0];
    for (size_t i = 0; i < v.size(); ++i) if (Pass(base, v[i])) mp.idx.push_back((int)i);
    const size_t n = mp.idx.size();
    mp.pp.resize(4 * n); mp.ee.resize(4 * n);
    auto add = [&](float* o, double p, double th, double ph, double m) {
        const double st = std::sin(th * d2r);
        o[0] += (float)(p * st * std::cos(ph * d2r)); o[1] += (float)(p * st * std::sin(ph * d2r));
        o[2] += (float)(p * std::cos(th * d2r)); o[3] += (float)std::sqrt(p * p + m * m);
    };
    for (size_t k = 0; k < n; ++k) {
        const Row& r = v[mp.idx[k]];
        float* a = &mp.pp[4 * k]; float* b = &mp.ee[4 * k];
        a[0] = a[1] = a[2] = a[3] = 0; b[0] = b[1] = b[2] = b[3] = 0;
        add(a, r.pip_p, r.pip_theta, r.pip_phi, mPi); add(a, r.pim_p, r.pim_theta, r.pim_phi, mPi);
        add(b, r.ep_p, r.ep_theta, r.ep_phi, mE); add(b, r.em_p, r.em_theta, r.em_phi, mE);
    }
}
inline void Kin(const MixPool& mp, size_t k, size_t j, MixKin& q, bool withEta) {
    const float* a = &mp.pp[4 * k]; const float* b = &mp.ee[4 * j];
    const double x = (double)a[0] + b[0], y = (double)a[1] + b[1], z = (double)a[2] + b[2], e = (double)a[3] + b[3];
    const double p2 = x * x + y * y + z * z, m2 = e * e - p2;
    q.m4 = 1e-3 * (m2 > 0 ? std::sqrt(m2) : -std::sqrt(-m2));
    const double pa = std::sqrt((double)a[0] * a[0] + (double)a[1] * a[1] + (double)a[2] * a[2]);
    const double pc = std::sqrt((double)b[0] * b[0] + (double)b[1] * b[1] + (double)b[2] * b[2]);
    double cs = ((double)a[0] * b[0] + (double)a[1] * b[1] + (double)a[2] * b[2]) / (pa * pc);
    q.oal = std::acos(std::max(-1., std::min(1., cs))) * TMath::RadToDeg();
    const double q2 = (mp.eIni - e) * (mp.eIni - e) - x * x - y * y - (mp.pb - z) * (mp.pb - z);
    q.mm = 1e-3 * (q2 > 0 ? std::sqrt(q2) : -std::sqrt(-q2));
    if (withEta) {
        const double mEta = 1e3 * scan::kMeta, eh = std::sqrt(p2 + mEta * mEta), g = eh / mEta;
        const double bx = x / eh, by = y / eh, bz = z / eh, b2 = bx * bx + by * by + bz * bz;
        double u[2][3];
        const float* w[2] = {a, b};
        for (int t = 0; t < 2; ++t) {
            const double bp = bx * w[t][0] + by * w[t][1] + bz * w[t][2];
            const double f = (b2 > 0 ? (g - 1) * bp / b2 : 0.) - g * w[t][3];
            u[t][0] = w[t][0] + f * bx; u[t][1] = w[t][1] + f * by; u[t][2] = w[t][2] + f * bz;
        }
        const double n0 = std::sqrt(u[0][0] * u[0][0] + u[0][1] * u[0][1] + u[0][2] * u[0][2]);
        const double n1 = std::sqrt(u[1][0] * u[1][0] + u[1][1] * u[1][1] + u[1][2] * u[1][2]);
        const double c2 = (u[0][0] * u[1][0] + u[0][1] * u[1][1] + u[0][2] * u[1][2]) / (n0 * n1);
        q.oae = std::acos(std::max(-1., std::min(1., c2))) * TMath::RadToDeg();
    }
}
// pass counts of selections on the mixed rows in [mlo, mhi): pas[half][sel][bin], tot[half][bin]
struct MixSpec { std::vector<std::vector<double>> pas[2]; std::vector<double> tot[2]; };
MixSpec MixedFill(const MixPool& mp, const std::vector<Sel>& sels, double mlo, double mhi, double* closure = nullptr) {
    MixSpec ms;
    for (int h = 0; h < 2; ++h) { ms.tot[h].assign(scan::kNbin, 0.); ms.pas[h].assign(sels.size(), std::vector<double>(scan::kNbin, 0.)); }
    const std::vector<Row>& v = *mp.rows;
    const size_t n = mp.idx.size();
    MixKin q;
    if (closure) {     // recomputation checked on unmixed rows
        double dm = 0, dl = 0, de = 0, dmm = 0;
        for (size_t k = 0; k < n; k += 7) {
            const Row& r = v[mp.idx[k]];
            Kin(mp, k, k, q, true);
            dm = std::max(dm, std::fabs(q.m4 - r.m4)); dl = std::max(dl, std::fabs(q.oal - r.oa_lab));
            de = std::max(de, std::fabs(q.oae - r.oa_eta)); if (r.mm4 > 0.5) dmm = std::max(dmm, std::fabs(q.mm - r.mm4));
        }
        closure[0] = dm; closure[1] = dl; closure[2] = de; closure[3] = dmm;
    }
    for (size_t k = 0; k < n; ++k) {
        const Row& a = v[mp.idx[k]];
        for (int o = 0; o < kNOff; ++o) {
            const size_t j = (k + kOff[o]) % n;
            const Row& b = v[mp.idx[j]];
            if (a.grp == b.grp) continue;
            Kin(mp, k, j, q, false);
            if (!(q.m4 >= mlo && q.m4 < mhi)) continue;
            Kin(mp, k, j, q, true);
            const int ib = (int)((q.m4 - scan::kMlo) / scan::kBinW);
            if (ib < 0 || ib >= scan::kNbin) continue;
            Row m = a;
            m.oa = b.oa; m.mee = b.mee; m.m4 = (float)q.m4; m.oa_lab = (float)q.oal; m.mm4 = (float)q.mm; m.oa_eta = (float)q.oae;
            m.ep_isUsedVertex = b.ep_isUsedVertex; m.em_isUsedVertex = b.em_isUsedVertex;
            const int h = o & 1;
            ms.tot[h][ib] += 1;
            for (size_t s = 0; s < sels.size(); ++s) if (sels[s](m)) ms.pas[h][s][ib] += 1;
        }
    }
    return ms;
}
std::vector<double> AddV(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> c(a.size()); for (size_t i = 0; i < a.size(); ++i) c[i] = a[i] + b[i]; return c;
}

// ---- turn-on interval (copy of scan_peak.C TurnOn) ----------------------------------
struct Turn { bool ok = false; double fmin = 0, fmax = 0, delta = 0, mA = -1, mB = -1, fLo = 0, fHi = 0, dfit = 0; bool edgeIn = false; };
Turn TurnOn(const std::vector<double>& pas, const std::vector<double>& tot, int reb, double minRows) {
    Turn t;
    std::vector<double> x, f;
    for (int i = 0; i + reb <= scan::kNbin; i += reb) {
        const double xc = scan::kMlo + (i + 0.5 * reb) * scan::kBinW;
        if (xc < 0.40 || xc > 0.80) continue;
        double a = 0, b = 0;
        for (int j = 0; j < reb; ++j) { a += pas[i + j]; b += tot[i + j]; }
        if (b < minRows) continue;
        x.push_back(xc); f.push_back(a / b);
    }
    const int n = (int)x.size();
    if (n < 5) return t;
    t.ok = true;
    t.fmin = *std::min_element(f.begin(), f.end()); t.fmax = *std::max_element(f.begin(), f.end());
    t.delta = t.fmax - t.fmin;
    double fl = 1e9, fh = -1e9; bool first = true;
    for (int i = 0; i < n; ++i) {
        if (x[i] < scan::kFitLo || x[i] > scan::kFitHi) continue;
        if (first) { t.fLo = f[i]; first = false; }
        t.fHi = f[i];
        fl = std::min(fl, f[i]); fh = std::max(fh, f[i]);
    }
    t.dfit = fh > fl ? fh - fl : 0;
    if (t.delta < 1e-6) return t;
    double mA = 1e9, mB = -1e9;
    for (int i = 0; i + 1 < n; ++i)
        for (double L : {0.1, 0.9}) {
            const double g0 = (f[i] - t.fmin) / t.delta - L, g1 = (f[i + 1] - t.fmin) / t.delta - L;
            if (g0 * g1 > 0 || g0 == g1) continue;
            const double m = x[i] + (x[i + 1] - x[i]) * g0 / (g0 - g1);
            mA = std::min(mA, m); mB = std::max(mB, m);
        }
    if (mA < 1e8) { t.mA = mA; t.mB = mB; }
    const bool aIn = t.mA >= scan::kFitLo && t.mA <= scan::kFitHi, bIn = t.mB >= scan::kFitLo && t.mB <= scan::kFitHi;
    t.edgeIn = (aIn || bIn) && t.dfit >= 0.05;
    return t;
}

// ---- fake peak (copy of scan_peak.C GaussCoef / FakeOne / FakeS) ---------------------
double GaussCoef(const std::vector<double>& y, const std::vector<double>& w, std::vector<double>* coef = nullptr) {
    std::vector<double> F, yy, var;
    const double x0 = 0.5 * (scan::kFitLo + scan::kFitHi), dx = 0.5 * (scan::kFitHi - scan::kFitLo);
    int npts = 0;
    for (int i = 0; i < scan::kNbin; ++i) {
        const double x1 = scan::kMlo + i * scan::kBinW, xc = x1 + 0.5 * scan::kBinW;
        if (xc < scan::kFitLo || xc >= scan::kFitHi) continue;
        const double t = (xc - x0) / dx;
        F.push_back(scan::GaussFrac(x1, x1 + scan::kBinW, scan::kPeakMean, scan::kPeakSigma));
        F.push_back(1.); F.push_back(t); F.push_back(t * t);
        yy.push_back(y[i]); var.push_back(std::max(w[i], 1e-9));
        npts++;
    }
    const scan::LinRes r = scan::LinFit(npts, 4, yy.data(), var.data(), F.data());
    if (!r.ok) return 0;
    if (coef) {
        coef->assign(scan::kNbin, 0.);
        int q = 0;
        for (int i = 0; i < scan::kNbin; ++i) {
            const double xc = scan::kMlo + (i + 0.5) * scan::kBinW;
            if (xc < scan::kFitLo || xc >= scan::kFitHi) continue;
            double a = 0;
            for (int l = 0; l < 4; ++l) a += r.cov[0][l] * F[q * 4 + l] / var[q];
            (*coef)[i] = a; q++;
        }
    }
    return r.par[0];
}
double FakeOne(const PeakFit& smooth, const std::vector<double>& pas, const std::vector<double>& tot, double* ebin = nullptr) {
    std::vector<double> y(scan::kNbin, 0.), w(scan::kNbin, 1.), ef(scan::kNbin, 0.), coef;
    for (int i = 0; i < scan::kNbin; ++i) {
        const double xc = scan::kMlo + (i + 0.5) * scan::kBinW;
        if (xc < scan::kFitLo || xc >= scan::kFitHi) continue;
        const double b = smooth.Bkg(xc);
        const double f = tot[i] > 0 ? pas[i] / tot[i] : 0;
        y[i] = b * f; w[i] = std::max(y[i], 1.);
        ef[i] = tot[i] > 0 ? b * std::sqrt(std::max(f * (1 - f), 0.) / tot[i]) : 0;
    }
    const double s = GaussCoef(y, w, &coef);
    if (ebin) { double e2 = 0; for (int i = 0; i < scan::kNbin; ++i) e2 += coef[i] * coef[i] * ef[i] * ef[i]; *ebin = std::sqrt(e2); }
    return s;
}
struct Fake { double s = 0, e = 0; };
Fake FakeS(const PeakFit& smooth, const std::vector<double>& p0, const std::vector<double>& t0,
           const std::vector<double>& p1, const std::vector<double>& t1) {
    Fake f; double eb = 0;
    f.s = FakeOne(smooth, AddV(p0, p1), AddV(t0, t1), &eb);
    const double a = FakeOne(smooth, p0, t0), b = FakeOne(smooth, p1, t1);
    f.e = std::max(eb, 0.5 * std::fabs(a - b));
    return f;
}

} // namespace optsc

// ============================================================================
void opt_sculpt() {
    using namespace opt;
    using namespace optsc;
    TStopwatch sw;
    gSystem->mkdir(kFigOut, kTRUE); gSystem->mkdir(kTabOut, kTRUE);
    Style();
    scan::Data d;
    if (!scan::Load(d, scan::kCompact, scan::kMlo, scan::kMhi)) return;
    const Cuts base = Pre();

    // ---- the tests ----
    std::vector<Test> T;
    auto add = [&](const char* fam, const std::string& lab, double x, const Cuts& c) { Test t; t.fam = fam; t.lab = lab; t.c = c; t.x = x; T.push_back(t); };
    for (int i = 1; i <= 8; ++i) { Cuts c = base; c.oaLo = 0.5f * i; add("OA >", Form("%.1f", 0.5 * i), 0.5 * i, c); }
    for (double v : {4., 6., 8., 10., 15., 20., 30.}) { Cuts c = base; c.oaHi = v; add("OA <", Form("%.0f", v), v, c); }
    for (int v = 25; v <= 70; v += 5) { Cuts c = base; c.oaLab = v; add("oa_lab <", Form("%d", v), v, c); }
    for (double v : {0.28, 0.30, 0.32, 0.34, 0.36, 0.38}) { Cuts c = base; c.mpLo = v; add("M(pipi) >", Form("%.2f", v), v, c); }
    for (double v : {0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50}) { Cuts c = base; c.mpHi = v; add("M(pipi) <", Form("%.2f", v), v, c); }
    for (int i = 0; i <= 7; ++i) { const double v = 2.60 + 0.05 * i; Cuts c = base; c.mm4Hi = v; add("MM4 <", Form("%.2f", v), v, c); }
    for (double v : {100., 120., 140., 150., 160.}) { Cuts c = base; c.etaLo = v; add("eta_rest >", Form("%.0f", v), v, c); }
    for (double v : {10., 15., 20., 30., 40., 60., 100.}) { Cuts c = base; c.dzPi = v; add("dz_pi <", Form("%.0f", v), v, c); }
    for (int v : {1, 2, 4}) { Cuts c = base; c.vtx = v; add("vertex fit", v == 1 ? ">=1 pi" : (v == 2 ? "both pi" : "all 4"), v, c); }
    add("named", "peak nominal", 0, PeakNom());
    add("named", "peak nominal + vertex", 0, PeakNomV());
    add("named", "DIRECT (frozen)", 0, Direct());
    add("named", "PEAK (frozen)", 0, PeakF());
    add("named", "pre + vertex", 0, PreV());
    add("named", "bump nominal", 0, BumpNom());
    add("named", "bump nominal + vertex", 0, BumpNomV());
    add("named", "slides", 0, Slides());
    add("named", "old code", 0, OldCode());
    add("named", "X17 band (frozen)", 0, X17());
    const int nt = (int)T.size();

    // ---- unlike-sign and like-sign spectra: pass and base ----
    std::vector<Sel> sels;
    sels.push_back(MkSel(base));
    for (const Test& t : T) sels.push_back(MkSel(t.c));
    const std::vector<Spectra> sp = scan::FillMany(d, sels);
    const PeakFit smooth = scan::FitPeakPois(sp[0]);      // Bref = its pol2 (CB inside)
    const PeakFit pm0 = scan::FitPeakMain(sp[0]);
    printf("spectra at %.1f s; base S = %.0f +- %.0f\n", sw.RealTime(), pm0.S, pm0.eS); sw.Continue();

    // ---- mixed rows ----
    MixPool mp; BuildPool(d, base, mp);
    std::vector<Sel> msels;
    for (const Test& t : T) msels.push_back(MkSel(t.c));
    double clos[4];
    const MixSpec mfit = MixedFill(mp, msels, scan::kFitLo, scan::kFitHi, clos);
    printf("mixed fit range at %.1f s (pool %zu)\n", sw.RealTime(), mp.idx.size()); sw.Continue();
    // full range for the pass-fraction figure: four cuts
    std::vector<int> figIdx;
    for (int i = 0; i < nt; ++i)
        if ((T[i].fam == "M(pipi) <" && T[i].lab == "0.42") || (T[i].fam == "MM4 <" && T[i].lab == "2.85") ||
            (T[i].fam == "eta_rest >" && T[i].lab == "140") || (T[i].fam == "oa_lab <" && T[i].lab == "55") ||
            (T[i].fam == "dz_pi <" && T[i].lab == "40") || (T[i].fam == "vertex fit" && T[i].lab == "both pi")) figIdx.push_back(i);
    std::vector<Sel> fsels;
    for (int i : figIdx) fsels.push_back(MkSel(T[i].c));
    const MixSpec mfull = MixedFill(mp, fsels, 0.40, 0.80);
    printf("mixed full range at %.1f s\n", sw.RealTime()); sw.Continue();

    // ---- evaluation ----
    struct Res { Turn ls, us, mx; Fake fk; PeakFit lsPP, lsMM; bool G = false, M = false; };
    std::vector<Res> R(nt);
    const double sigS = pm0.eS;
    const std::vector<double> lsTot = AddV(sp[0].n[1], sp[0].n[2]);
    for (int i = 0; i < nt; ++i) {
        const Spectra& s = sp[1 + i];
        R[i].ls = TurnOn(AddV(s.n[1], s.n[2]), lsTot, 8, 20);
        R[i].us = TurnOn(s.n[0], sp[0].n[0], 4, 50);
        R[i].mx = TurnOn(AddV(mfit.pas[0][i], mfit.pas[1][i]), AddV(mfit.tot[0], mfit.tot[1]), 4, 50);
        R[i].fk = FakeS(smooth, mfit.pas[0][i], mfit.tot[0], mfit.pas[1][i], mfit.tot[1]);
        R[i].lsPP = scan::FitPeakLikeSign(s, 1); R[i].lsMM = scan::FitPeakLikeSign(s, 2);
        R[i].G = R[i].ls.edgeIn;
        R[i].M = std::fabs(R[i].fk.s) > 0.25 * sigS;
    }
    // method check: the mixed spectrum itself, fitted with Gaussian + pol2, at the data normalisation
    double mixSelf = 0, mixSelfE = 0;
    {
        std::vector<double> y(scan::kNbin, 0.), w(scan::kNbin, 1.), coef;
        const std::vector<double> tot = AddV(mfit.tot[0], mfit.tot[1]);
        double nmix = 0, ndat = 0;
        for (int i = 0; i < scan::kNbin; ++i) {
            const double xc = scan::kMlo + (i + 0.5) * scan::kBinW;
            if (xc < scan::kFitLo || xc >= scan::kFitHi) continue;
            y[i] = tot[i]; w[i] = std::max(tot[i], 1.); nmix += tot[i]; ndat += smooth.Bkg(xc);
        }
        const double s = GaussCoef(y, w, &coef);
        double e2 = 0; for (int i = 0; i < scan::kNbin; ++i) e2 += coef[i] * coef[i] * w[i];
        mixSelf = s * ndat / nmix; mixSelfE = std::sqrt(e2) * ndat / nmix;
    }

    // ---- table ----
    Out o;
    o.Open(std::string(kTabOut) + "/opt_sculpt.txt");
    o.P("# opt_sculpt.txt  --  written by plots/exp_report/opt_sculpt.C\n");
    o.P("# EXPERIMENTAL DATA. Base = pre: cor_ok == 1, isBest >= 0, vz > -160 mm (FAT: vz > -500, PT3). Corrected momenta.\n");
    o.P("# Sculpting test of every cut value: pass fraction f against M4 for US (unlike-sign rows), LS (like-sign, ++ and -- added),\n");
    o.P("# MIX (pion pair of one row + lepton pair of another event group; 20 offsets; pool = %zu base rows).\n", mp.idx.size());
    o.P("# Recomputation check on unmixed rows: max |dM4| = %.1e GeV, |d oa_lab| = %.1e deg, |d eta_rest| = %.1e deg, |d MM4| = %.1e GeV.\n",
        clos[0], clos[1], clos[2], clos[3]);
    o.P("# f lo / f hi = f at the two ends of the fit range %.4f-%.4f, dfit = change of f inside it; M10-M90 = where f crosses 10 %% / 90 %% of its\n", scan::kFitLo, scan::kFitHi);
    o.P("# total change in 0.40-0.80 (LS 20 MeV bins, US 10 MeV bins; MIX filled in the fit range only: M10-M90 not computed, printed as -1).\n");
    o.P("# S_fake = Gaussian content of Bref x f_MIX fitted with the frozen Gaussian + pol2 (Bref = pol2 of the Poisson fit of pre, CB inside),\n");
    o.P("#   in rows at the data normalisation; err = max(binomial, half difference of the two offset halves).\n");
    o.P("# Flags: G = LS turn-on inside the fit range with dfit >= 0.05; M = |S_fake| > 0.25 sigma_S(pre) = %.0f rows. G or M: the value can bias S.\n", 0.25 * sigS);
    o.P("# Like-sign null test: Poisson fit of N++ and N-- with the frozen Gaussian + pol2 (expected 0).\n");
    o.P("# Base: S(pre) = %.0f +- %.0f (P_main); method check: the mixed spectrum itself fitted with Gaussian + pol2 = %.0f +- %.0f rows (data normalisation).\n#\n",
        pm0.S, pm0.eS, mixSelf, mixSelfE);
    o.P("%-11s %-22s | %5s %5s %5s %11s | %5s %5s %5s %11s | %5s %5s %5s %11s | %6s %4s | %6s %4s %6s %4s | G M\n", "cut", "value",
        "LS lo", "f hi", "dfit", "M10-M90", "US lo", "f hi", "dfit", "M10-M90", "MX lo", "f hi", "dfit", "M10-M90", "S_fake", "err", "S++", "err", "S--", "err");
    std::string last;
    for (int i = 0; i < nt; ++i) {
        if (T[i].fam != last && !last.empty()) o.P("\n");
        last = T[i].fam;
        const Res& r = R[i];
        o.P("%-11s %-22s | %5.3f %5.3f %5.3f %5.3f-%5.3f | %5.3f %5.3f %5.3f %5.3f-%5.3f | %5.3f %5.3f %5.3f %5.3f-%5.3f | %6.0f %4.0f | %6.0f %4.0f %6.0f %4.0f | %d %d\n",
            T[i].fam.c_str(), T[i].lab.c_str(), r.ls.fLo, r.ls.fHi, r.ls.dfit, r.ls.mA, r.ls.mB, r.us.fLo, r.us.fHi, r.us.dfit, r.us.mA, r.us.mB,
            r.mx.fLo, r.mx.fHi, r.mx.dfit, -1.0, -1.0, r.fk.s, r.fk.e, r.lsPP.S, r.lsPP.eS, r.lsMM.S, r.lsMM.eS, (int)r.G, (int)r.M);
    }
    o.P("\nrun time %.1f s\n", sw.RealTime());
    o.Close();

    // ---- figure 1: pass fractions against M4 ----
    {
        TCanvas* c = new TCanvas("c_pf", "", 1500, 1000);
        c->Divide(3, 2, 0.001, 0.001);
        const char* nm[6] = {"oa_lab < 55#circ", "M(#pi^{+}#pi^{-}) < 0.42 GeV/c^{2}", "MM_{4} < 2.85 GeV/c^{2} (frozen PEAK)",
                             "#eta-rest angle > 140#circ", "#Deltaz_{#pi} < 40 mm (frozen)", "both #pi in the vertex fit"};
        const char* fam[6] = {"oa_lab <", "M(pipi) <", "MM4 <", "eta_rest >", "dz_pi <", "vertex fit"};
        for (int q = 0; q < 6; ++q) {
            int i = -1, k = -1;
            for (size_t m = 0; m < figIdx.size(); ++m) if (T[figIdx[m]].fam == fam[q]) { i = figIdx[m]; k = (int)m; }
            if (i < 0) continue;
            TPad* p = (TPad*)c->cd(q + 1);
            p->SetLeftMargin(0.16); p->SetRightMargin(0.03); p->SetTopMargin(0.09); p->SetBottomMargin(0.15);
            TH1D* fr = Frame(Form("fpf%d", q), 0.40, 0.80, 0, 1.15, "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]", "pass fraction");
            fr->GetYaxis()->SetNdivisions(506);
            TBox* b = new TBox(scan::kFitLo, 0, scan::kFitHi, 1.15); b->SetFillColor(TColor::GetColor("#ececec")); b->SetLineWidth(0); b->Draw();
            TBox* bw = new TBox(scan::kPeakLo, 0, scan::kPeakHi, 1.15); bw->SetFillColor(TColor::GetColor("#d9ddf6")); bw->SetLineWidth(0); bw->Draw();
            const Spectra& s = sp[1 + i];
            TGraphErrors* gu = GE(kBlack, 20, 0.6); TGraphErrors* gl = GE(kRed + 1, 21, 0.7); TGraphErrors* gm = GE(kBlue + 1, 1, 0.1);
            gm->SetLineWidth(3);
            for (int j = 0; j + 4 <= scan::kNbin; j += 4) {
                double a = 0, t = 0, al = 0, tl = 0, am = 0, tm = 0;
                for (int u = 0; u < 4; ++u) {
                    a += s.n[0][j + u]; t += sp[0].n[0][j + u];
                    am += mfull.pas[0][k][j + u] + mfull.pas[1][k][j + u]; tm += mfull.tot[0][j + u] + mfull.tot[1][j + u];
                }
                for (int u = 0; u < 8 && j + u < scan::kNbin; ++u) { if (j % 8) break; al += s.n[1][j + u] + s.n[2][j + u]; tl += lsTot[j + u]; }
                const double xc = scan::kMlo + (j + 2) * scan::kBinW;
                if (xc < 0.40 || xc > 0.80) continue;
                if (t > 50) AddPt(gu, xc, a / t, std::sqrt(std::max(a / t * (1 - a / t), 0.) / t));
                if (tm > 50) AddPt(gm, xc, am / tm, 0);
                if (j % 8 == 0 && tl > 20) AddPt(gl, xc + 2 * scan::kBinW, al / tl, std::sqrt(std::max(al / tl * (1 - al / tl), 0.) / tl));
            }
            gm->Draw("L SAME"); gl->Draw("PZ SAME"); gu->Draw("PZ SAME");
            Txt(0.16, 0.93, Form("(%c) %s", 'a' + q, nm[q]), 0.055);
            Txt(0.20, 0.20, Form("S_{fake} = %.0f #pm %.0f rows", R[i].fk.s, R[i].fk.e), 0.050, R[i].M ? kRed + 1 : kBlue + 2);
            if (q == 0) {
                TLegend* l = Leg(0.45, 0.29, 0.97, 0.52, 0.045);
                l->AddEntry(gu, "e^{+}e^{-} rows (US)", "pe"); l->AddEntry(gl, "e^{#pm}e^{#pm} rows (LS)", "pe");
                l->AddEntry(gm, "mixed rows (MIX)", "l"); l->Draw();
            }
            if (q == 1) Txt(0.20, 0.27, "grey: fit range; blue: peak window", 0.042, kGray + 2);
        }
        Save(c, "opt_sculpt_passfrac");
    }
    // ---- figure 2: S_fake against the cut value ----
    {
        TCanvas* c = new TCanvas("c_fk", "", 1500, 1000);
        c->Divide(3, 2, 0.001, 0.001);
        const char* fam[6] = {"oa_lab <", "M(pipi) >", "M(pipi) <", "MM4 <", "eta_rest >", "dz_pi <"};
        const char* ax[6] = {"upper cut on oa_{lab}((#pi#pi),(ee)) [deg]", "lower cut on M(#pi^{+}#pi^{-}) [GeV/c^{2}]",
                             "upper cut on M(#pi^{+}#pi^{-}) [GeV/c^{2}]", "upper cut on MM_{4} [GeV/c^{2}]",
                             "lower cut on the #eta-rest angle [deg]", "upper cut on #Deltaz_{#pi} [mm]"};
        for (int q = 0; q < 6; ++q) {
            TPad* p = (TPad*)c->cd(q + 1);
            p->SetLeftMargin(0.17); p->SetRightMargin(0.03); p->SetTopMargin(0.09); p->SetBottomMargin(0.15);
            double x1 = 1e9, x2 = -1e9;
            TGraphErrors* g = GE(kBlue + 1, 20, 1.1); TGraphErrors* gf = GE(kRed + 1, 24, 1.4);
            for (int i = 0; i < nt; ++i) {
                if (T[i].fam != fam[q]) continue;
                x1 = std::min(x1, T[i].x); x2 = std::max(x2, T[i].x);
                AddPt(g, T[i].x, R[i].fk.s, R[i].fk.e);
                if (R[i].G || R[i].M) AddPt(gf, T[i].x, R[i].fk.s, 0);
            }
            const double dx = 0.08 * (x2 - x1);
            TH1D* fr = Frame(Form("ffk%d", q), x1 - dx, x2 + dx, -500, 500, ax[q], "S_{fake} [rows]");
            fr->GetYaxis()->SetNdivisions(505);
            TBox* b = new TBox(x1 - dx, -0.25 * sigS, x2 + dx, 0.25 * sigS); b->SetFillColor(TColor::GetColor("#cdf2cd")); b->SetLineWidth(0); b->Draw();
            HLine(x1 - dx, x2 + dx, 0, kGray + 2, 2, 1);
            g->Draw("PZ SAME"); if (gf->GetN()) gf->Draw("P SAME");
            fr->Draw("AXIS SAME");
            const char* nice[6] = {"oa_{lab} < X", "M(#pi^{+}#pi^{-}) > X", "M(#pi^{+}#pi^{-}) < X", "MM_{4} < X", "#eta-rest angle > X", "#Deltaz_{#pi} < X"};
            Txt(0.17, 0.93, Form("(%c) %s", 'a' + q, nice[q]), 0.055);
            // points outside the frame
            for (int i = 0; i < nt; ++i) {
                if (T[i].fam != fam[q] || std::fabs(R[i].fk.s) < 500) continue;
                TLatex* t = new TLatex(T[i].x, R[i].fk.s > 0 ? 420 : -420, Form("#%s %.0f", R[i].fk.s > 0 ? "uparrow" : "downarrow", R[i].fk.s));
                t->SetTextAlign(R[i].fk.s > 0 ? 33 : 31); t->SetTextSize(0.045); t->SetTextFont(42); t->SetTextColor(kRed + 1); t->Draw();
            }
            if (q == 0) {
                TLegend* l = Leg(0.22, 0.17, 0.95, 0.36, 0.045);
                l->AddEntry(g, "S_{fake} from the mixed rows", "pe");
                l->AddEntry(gf, "flagged (G or M)", "p");
                l->AddEntry(b, Form("#pm0.25 #sigma_{S}(pre) = #pm%.0f rows", 0.25 * sigS), "f");
                l->Draw();
            }
        }
        Save(c, "opt_sculpt_fake");
    }
    printf("opt_sculpt: done in %.1f s\n", sw.RealTime());
}
