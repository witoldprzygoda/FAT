// ==============================================================================================================
// vtx_lib.h -- copy (2026-11-14, integration of the event-vertex efficiency into report_all / report_exp) of
//   /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/vtx_eff/ve_verify/vv_lib.h
//   (study of the eta efficiency of the mandatory event-vertex cut vz > -500 mm; verifier ve_verify, own code).
// What it does: own spectra, like-sign CB, GLS peak fit (P_main convention), edge-step fit and eps = X_V/(X_V+X_Nv).
// Changes against the original: only the file name and the include names (vv_* -> vtx_*); the code is unchanged.
// The original outputs are copied to docs/latex/tables_all/vtx_*.txt; the notes to docs/latex/notes_all/ve_*_note.md.
// Run (from /home/przygoda/HADES/FAT/FAT; vtx_fill.C first, then the analysis macros, at most 6 processes at once):
//   B='gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/vtx_eff/ve_verify/aclic", kTRUE)'
//   (header; included by the vtx_*.C macros below)
// ==============================================================================================================
// vtx_lib.h -- ve_verify: own spectra, CB and peak/edge fits (written from scratch; same conventions as the report:
// 2.5 MeV fine bins 0.35-0.95, CB = 2 sqrt(N++ N--) per 10 MeV linearly interpolated to the fine bins,
// P_main = GLS fit of N - CB with frozen Gaussian (0.5439, 5.4 MeV) + pol2 in 0.5225-0.61, covariance = expected
// Poisson (iterated) + covariance of the interpolated CB).
#ifndef VV_LIB_H
#define VV_LIB_H
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>
#include <algorithm>
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TDecompChol.h"
#include "TMath.h"

namespace vvl {

struct R {             // one row in memory
  int file, grp, smp, day, best, start, vcand, ng;
  float vz, m4, mee, mpp, oa, meeg, m5c, zx, rx, dca, zp, zxe, rxe;
  signed char nwin, cor;
};

const int kNb = 240; const double kLo = 0.35, kHi = 0.95, kBw = 0.0025; const int kReb = 4;

struct Spec {
  std::vector<double> n[3], w[3], cbW, cbWv, cb;
  double winN = 0, winVarG = 0;        // us rows in the peak window and their event-group variance
  long long nsel[3] = {0, 0, 0};
  int nw = kNb / kReb;
  double x(int i) const { return kLo + (i + 0.5) * kBw; }
  void Interp(int i, int& j0, int& j1, double& f) const {
    const double t = (x(i) - (kLo + 0.5 * kReb * kBw)) / (kReb * kBw);
    if (t <= 0) { j0 = j1 = 0; f = 0; return; }
    if (t >= nw - 1) { j0 = j1 = nw - 1; f = 0; return; }
    j0 = (int)t; j1 = j0 + 1; f = t - j0;
  }
  double CbCov(int i, int k) const {
    int a0, a1, b0, b1; double fa, fb; Interp(i, a0, a1, fa); Interp(k, b0, b1, fb);
    const int ja[2] = {a0, a1}, jb[2] = {b0, b1};
    const double wa[2] = {(1 - fa) / kReb, fa / kReb}, wb[2] = {(1 - fb) / kReb, fb / kReb};
    double s = 0;
    for (int p = 0; p < 2; ++p) for (int q = 0; q < 2; ++q) if (ja[p] == jb[q]) s += wa[p] * wb[q] * cbWv[ja[p]];
    return s;
  }
  void Finish(double k = 1.0) {
    cbW.assign(nw, 0); cbWv.assign(nw, 0); cb.assign(kNb, 0);
    for (int s = 0; s < 3; ++s) { w[s].assign(nw, 0); for (int i = 0; i < kNb; ++i) w[s][i / kReb] += n[s][i]; }
    for (int j = 0; j < nw; ++j) { cbW[j] = k * 2 * std::sqrt(w[1][j] * w[2][j]); cbWv[j] = k * k * (w[1][j] + w[2][j]); }
    for (int i = 0; i < kNb; ++i) { int j0, j1; double f; Interp(i, j0, j1, f); cb[i] = ((1 - f) * cbW[j0] + f * cbW[j1]) / kReb; }
  }
};

typedef std::function<bool(const R&)> Sel;
typedef std::function<float(const R&)> Var;

// fill several selections in one pass; v = histogrammed variable (default m4); peak window for the group variance
inline std::vector<Spec> Fill(const std::vector<R>& rows, const std::vector<Sel>& sels, Var v = nullptr,
                              double wlo = 0.53364, double whi = 0.55416) {
  const int ns = sels.size();
  std::vector<Spec> out(ns);
  for (auto& s : out) for (int k = 0; k < 3; ++k) s.n[k].assign(kNb, 0);
  std::vector<int> cnt(ns, 0);
  int cf = -1, cg = -1;
  for (const R& r : rows) {
    if (r.file != cf || r.grp != cg) {
      for (int s = 0; s < ns; ++s) if (cnt[s]) { out[s].winVarG += (double)cnt[s] * cnt[s]; cnt[s] = 0; }
      cf = r.file; cg = r.grp;
    }
    if (!r.cor) continue;
    const double m = v ? v(r) : r.m4;
    if (!(m == m)) continue;
    const int ib = (m >= kLo && m < kHi) ? (int)((m - kLo) / kBw) : -1;
    for (int s = 0; s < ns; ++s) {
      if (!sels[s](r)) continue;
      out[s].nsel[r.smp]++;
      if (ib >= 0) out[s].n[r.smp][ib] += 1;
      if (r.smp == 0 && m >= wlo && m < whi) { out[s].winN += 1; cnt[s]++; }
    }
  }
  for (int s = 0; s < ns; ++s) if (cnt[s]) out[s].winVarG += (double)cnt[s] * cnt[s];
  for (auto& s : out) s.Finish();
  return out;
}

inline double GFrac(double a, double b, double mu, double sg) {
  const double q = 1 / (std::sqrt(2.) * sg); return 0.5 * (std::erf((b - mu) * q) - std::erf((a - mu) * q));
}

// ---- generic GLS fit of y = F p on fine bins idx; var = max(model + add, 1) iterated; + CB covariance ----
struct Fit { bool ok = false; std::vector<double> p; TMatrixD cov; double chi2 = 0; int ndf = 0; double prob = 0;
             double S() const { return p.empty() ? 0 : p[0]; } double eS() const { return p.empty() ? 0 : std::sqrt(cov(0, 0)); } };

inline Fit GLS(const Spec& sp, const std::vector<int>& idx, const std::vector<std::vector<double>>& F,
               const std::vector<double>& y, const std::vector<double>& add, bool cbcov, double cbScale2 = 1.0, int sIdx = 0) {
  Fit r; const int n = idx.size(), np = F[0].size();
  if (n <= np) return r;
  TMatrixD C0(n, n);
  if (cbcov) for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) C0(a, b) = cbScale2 * sp.CbCov(idx[a], idx[b]);
  std::vector<double> var(n);
  for (int a = 0; a < n; ++a) var[a] = std::max(y[a] + add[a], 1.0);
  TMatrixD A(n, np); TVectorD Y(n);
  for (int a = 0; a < n; ++a) { Y[a] = y[a]; for (int k = 0; k < np; ++k) A(a, k) = F[a][k]; }
  TVectorD p(np); TMatrixD covp(np, np); double chi2 = 0;
  for (int it = 0; it < 8; ++it) {
    TMatrixD C(C0); for (int a = 0; a < n; ++a) C(a, a) += var[a];
    TMatrixD Ci(C); double det; Ci.Invert(&det);
    TMatrixD At(TMatrixD::kTransposed, A);
    TMatrixD AtCi = At * Ci;
    TMatrixD M = AtCi * A; M.Invert(&det);
    TVectorD b = AtCi * Y;
    TVectorD pn = M * b;
    const double old = it ? p[sIdx] : 1e30;
    p = pn; covp = M;
    TVectorD res = Y - A * p; chi2 = res * (Ci * res);
    for (int a = 0; a < n; ++a) { double f = 0; for (int k = 0; k < np; ++k) f += A(a, k) * p[k]; var[a] = std::max(f + add[a], 1.0); }
    if (std::fabs(p[sIdx] - old) < 1e-4 * std::sqrt(covp(sIdx, sIdx))) break;
  }
  r.ok = true; r.p.assign(np, 0); for (int k = 0; k < np; ++k) r.p[k] = p[k];
  r.cov.ResizeTo(np, np); r.cov = covp; r.chi2 = chi2; r.ndf = n - np; r.prob = TMath::Prob(chi2, r.ndf);
  return r;
}

struct PCfg { double mu = 0.5439, sg = 0.0054, lo = 0.5225, hi = 0.61; int ord = 2; int what = 0; // 0: N-CB; 1: N++; 2: N--; 3: N++ + N--; 4: N with free CB scale
              double k = 1.0; };

// peak fit: S (whole Gaussian) + polynomial
inline Fit Peak(const Spec& sp0, const PCfg& c) {
  Spec sp = sp0; if (c.k != 1.0) sp.Finish(c.k);
  std::vector<int> idx;
  for (int i = 0; i < kNb; ++i) if (sp.x(i) >= c.lo && sp.x(i) < c.hi) idx.push_back(i);
  const double x0 = 0.5 * (c.lo + c.hi), dx = 0.5 * (c.hi - c.lo);
  std::vector<std::vector<double>> F; std::vector<double> y, add;
  for (int i : idx) {
    std::vector<double> row; row.push_back(GFrac(kLo + i * kBw, kLo + (i + 1) * kBw, c.mu, c.sg));
    const double t = (sp.x(i) - x0) / dx; double tp = 1;
    for (int j = 0; j <= c.ord; ++j) { row.push_back(tp); tp *= t; }
    if (c.what == 4) row.push_back(sp.cb[i]);
    F.push_back(row);
    if (c.what == 0) { y.push_back(sp.n[0][i] - sp.cb[i]); add.push_back(sp.cb[i]); }
    else if (c.what == 1) { y.push_back(sp.n[1][i]); add.push_back(0); }
    else if (c.what == 2) { y.push_back(sp.n[2][i]); add.push_back(0); }
    else if (c.what == 3) { y.push_back(sp.n[1][i] + sp.n[2][i]); add.push_back(0); }
    else if (c.what == 6) { y.push_back(sp.cb[i]); add.push_back(0); }   // the CB itself (its "peak" content)
    else { y.push_back(sp.n[0][i]); add.push_back(0); }                  // 4, 5: N
  }
  return GLS(sp, idx, F, y, add, c.what == 0 || c.what == 4 || c.what == 6);
}

// free position and width: profile chi2 over a grid
struct Free { double mu = 0, sg = 0, S = 0, eS = 0, chi2 = 1e30, emu = 0, esg = 0; bool edgeHit = false; };
inline Free FreeFit(const Spec& sp, PCfg c) {
  Free b; std::vector<double> mus, sgs, ch;
  for (double mu = 0.5340; mu <= 0.55401; mu += 0.00025)
    for (double sg = 0.0030; sg <= 0.01201; sg += 0.00025) {
      c.mu = mu; c.sg = sg; Fit f = Peak(sp, c);
      if (!f.ok) continue;
      mus.push_back(mu); sgs.push_back(sg); ch.push_back(f.chi2);
      if (f.chi2 < b.chi2) { b.chi2 = f.chi2; b.mu = mu; b.sg = sg; b.S = f.S(); b.eS = f.eS(); }
    }
  // 1-sigma ranges (profile, dchi2 < 1)
  double mlo = 1e9, mhi = -1e9, slo = 1e9, shi = -1e9;
  for (size_t i = 0; i < ch.size(); ++i) if (ch[i] < b.chi2 + 1) {
    mlo = std::min(mlo, mus[i]); mhi = std::max(mhi, mus[i]); slo = std::min(slo, sgs[i]); shi = std::max(shi, sgs[i]); }
  b.emu = 0.5 * (mhi - mlo); b.esg = 0.5 * (shi - slo);
  b.edgeHit = (b.mu < 0.5341 || b.mu > 0.5539 || b.sg < 0.00301 || b.sg > 0.01199);
  return b;
}

// window sums (peak window) of N, CB
inline void WinSums(const Spec& sp, double a, double b, double& N, double& CB) {
  N = 0; CB = 0;
  for (int i = 0; i < kNb; ++i) {
    const double x1 = kLo + i * kBw, x2 = x1 + kBw, o = std::min(x2, b) - std::max(x1, a);
    if (o <= 0) continue; const double u = o / kBw; N += u * sp.n[0][i]; CB += u * sp.cb[i];
  }
}

// ---- bump edge step (report Eq. edge): N - CB = pol + [eta] + h * step(m), edge fixed ----
const double kMeta = 0.547862, kMpi0 = 0.1349768, kScale = 0.5439 / 0.547862, kScalePiPi = 0.9970;
inline double M4max(double mpp) {
  const double E0 = (kMeta * kMeta + kMpi0 * kMpi0 - mpp * mpp) / (2 * kMeta), p2 = E0 * E0 - kMpi0 * kMpi0;
  if (!(p2 >= 0)) return -1; const double m2 = kMeta * kMeta - kMeta * (E0 - std::sqrt(p2)); return m2 > 0 ? std::sqrt(m2) : -1;
}
inline double EdgeMeas(double mpp) { const double e = M4max(mpp / kScalePiPi); return e > 0 ? kScale * e : -1; }
struct ECfg { double mpLo, mpHi; double half = 0.045; double s = 0.0054; int ord = 2; bool useEdge = false; double edge = 0; };
inline double StepT(double m, const ECfg& c) {
  const double q = 1 / (std::sqrt(2.) * c.s);
  if (c.useEdge) return 0.5 * std::erfc((m - c.edge) * q);
  double s = 0; int u = 0;
  for (int i = 0; i < 21; ++i) { const double mp = c.mpLo + (i + 0.5) * (c.mpHi - c.mpLo) / 21; const double e = EdgeMeas(mp);
    if (e <= 0) continue; s += 0.5 * std::erfc((m - e) * q); ++u; }
  return u ? s / u : 0;
}
inline Fit Edge(const Spec& sp, const ECfg& c, bool& withEta) {
  const double e = c.useEdge ? c.edge : EdgeMeas(0.5 * (c.mpLo + c.mpHi));
  double lo = e - c.half, hi = e + c.half;
  const double mu = 0.5439, sg = 0.0054;
  if (hi <= mu - 3 * sg) withEta = false;
  else if (hi < mu - sg) { hi = mu - 3 * sg; withEta = false; }
  else { hi = std::max(hi, mu + 3 * sg); withEta = true; }
  lo = kLo + kBw * std::floor((lo - kLo) / kBw + 0.5); hi = kLo + kBw * std::floor((hi - kLo) / kBw + 0.5);
  std::vector<int> idx; for (int i = 0; i < kNb; ++i) if (sp.x(i) >= lo && sp.x(i) < hi) idx.push_back(i);
  const double x0 = 0.5 * (lo + hi), dx = 0.5 * (hi - lo);
  std::vector<std::vector<double>> F; std::vector<double> y, add;
  for (int i : idx) {
    std::vector<double> row; row.push_back(StepT(sp.x(i), c));
    const double t = (sp.x(i) - x0) / dx; double tp = 1;
    for (int j = 0; j <= c.ord; ++j) { row.push_back(tp); tp *= t; }
    if (withEta) row.push_back(GFrac(sp.x(i) - 0.5 * kBw, sp.x(i) + 0.5 * kBw, mu, sg));
    F.push_back(row); y.push_back(sp.n[0][i] - sp.cb[i]); add.push_back(sp.cb[i]);
  }
  return GLS(sp, idx, F, y, add, true);
}

inline double EpsV(double a, double ea, double b, double eb, double& err) {   // a/(a+b), independent errors
  const double s = a + b; err = std::sqrt(b * b * ea * ea + a * a * eb * eb) / (s * s); return a / s;
}

}  // namespace vvl
#endif
