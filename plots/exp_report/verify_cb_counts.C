// ============================================================================
// verify_cb_counts.C  --  independent check of the main numbers of the chapter
//                         "Combinatorial background" (prefix cb)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Written by the verifier of the cb chapter. It does NOT use the compact file,
// ScanEngine.h or cb_common.h: it reads the FAT outputs (tree scan_nt of
// output_pippim{epem,epep,emem}_exp.root) and recomputes, with its own code:
//   1. per selection and window: N, N++, N--, CB = 2 sqrt(N++ N--) in 10 MeV
//      bins of M4 = M(pi+pi-e+e-) (corrected momenta), linearly interpolated
//      between the 10 MeV bin centres and integrated exactly over the window,
//      CB/N, (N++ + N--)/CB and W = N - CB with event-group errors (event group
//      = adjacent rows with equal file_id, runnumber, event and vertex);
//      kept fractions of the vertex cuts;
//   2. the frozen eta fit (Gaussian 0.5439 GeV/c2, sigma 5.4 MeV, + pol2,
//      0.5225-0.6100) on N++ + N-- with the ROOT binned Poisson likelihood
//      (TH1::Fit options "L I") -- null test;
//   3. P_main-like chi2 fit (ROOT, option "I") of N - CB with bin errors
//      sqrt(N + var CB) (no CB covariance between bins, unlike ScanEngine),
//      and the change of S for CB x 0.97 and x 1.03;
//   4. the class-3 isBest test: join of the three trees by (file_id,
//      runnumber, event, vertex x, y, z); class n = n(e+) + n(e-) of the event
//      (distinct lepton tracks); kept fractions of isBest == 1 in the bump window.
// Selections (rows with cor_ok == 1; FAT already applied vertex z > -500 mm and
// PT3 = trigbit 8192, both re-checked here):
//   fat = isBest >= 0 ; pre = fat + vz > -160 mm ; vertex = pip_isUsedVertex &&
//   pim_isUsedVertex ; OA>2 = oa > 2 deg ; peak nominal = pre + oa_lab < 55 +
//   M(pipi) > 0.30 ; bump nominal = pre + oa_lab < 40 + 0.30 < M(pipi) < 0.41 ;
//   DIRECT = pre + dz_pi < 40 mm (dz_pi = max(|pip_z - vz|, |pim_z - vz|)) ;
//   PEAK = DIRECT + MM4 < 2.85 GeV/c2 ; "_b1" = the same with isBest == 1.
// Windows: data peak 0.53364-0.55416, data bump 0.4550-0.5225 GeV/c2.
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/output_pippim{epem,epep,emem}_exp.root
// Outputs: printout (redirect to a log), e.g.
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/cb_verify/verify_cb_counts.log
// Run (from /home/przygoda/HADES/FAT/FAT), about 1 min, 2 GB of memory:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_cb_verify", kTRUE)' 'plots/exp_report/verify_cb_counts.C+' > /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/cb_verify/verify_cb_counts.log 2>&1
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TMath.h"
#include "TTree.h"

namespace vcb {

const char* const kDir = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530";
const double kLo = 0.35, kHi = 0.95, kBw = 0.0025;
const int kNb = 240, kNw = 60;                      // 2.5 MeV fine bins, 10 MeV CB bins
const double kPkLo = 0.5439 - 1.9 * 0.0054, kPkHi = 0.5439 + 1.9 * 0.0054;
const double kBmLo = 0.4550, kBmHi = 0.5225;

struct R {
    float m4, mpipi, mee, oa, oalab, mm4, vx, vy, vz, pipz, pimz, pipv, pimv, isBest, cor, trig, file, run, ev;
    float epth, epph, emth, emph;
    int grp;
};

std::vector<R> Read(const char* tag) {
    std::vector<R> v;
    TFile* f = TFile::Open(Form("%s/output_pippim%s_exp.root", kDir, tag));
    TTree* t = (TTree*)f->Get("scan_nt");
    R r;
    t->SetBranchStatus("*", 0);
    const char* nm[] = {"m_pippimepem_cor", "m_pippim_cor", "m_ee_cor", "oa", "oa_pippim_epem_lab_cor", "mm_pippimepem_cor",
                        "eVertReco_x", "eVertReco_y", "eVertReco_z", "pip_z", "pim_z", "pip_isUsedVertex", "pim_isUsedVertex",
                        "isBest", "cor_ok", "trigbit", "file_id", "runnumber", "event", "ep_theta", "ep_phi", "em_theta", "em_phi"};
    float* ad[] = {&r.m4, &r.mpipi, &r.mee, &r.oa, &r.oalab, &r.mm4, &r.vx, &r.vy, &r.vz, &r.pipz, &r.pimz, &r.pipv, &r.pimv,
                   &r.isBest, &r.cor, &r.trig, &r.file, &r.run, &r.ev, &r.epth, &r.epph, &r.emth, &r.emph};
    for (int i = 0; i < 23; ++i) { t->SetBranchStatus(nm[i], 1); t->SetBranchAddress(nm[i], ad[i]); }
    const Long64_t n = t->GetEntries();
    v.reserve(n);
    int g = -1; float pf = -1, pr = -1, pe = -1, pz = 1e9;
    long badTrig = 0, badZ = 0;
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (r.file != pf || r.run != pr || r.ev != pe || r.vz != pz) { ++g; pf = r.file; pr = r.run; pe = r.ev; pz = r.vz; }
        r.grp = g;
        if (r.trig != 8192) ++badTrig;
        if (!(r.vz > -500)) ++badZ;
        v.push_back(r);
    }
    printf("%s: %zu rows, %d event groups, rows with trigbit != 8192: %ld, with vz <= -500: %ld\n", tag, v.size(), g + 1, badTrig, badZ);
    f->Close();
    return v;
}

typedef std::function<bool(const R&)> Sel;
inline bool Fat(const R& r) { return r.cor == 1 && r.isBest >= 0; }
inline bool Pre(const R& r) { return Fat(r) && r.vz > -160.f; }
inline bool Vtx(const R& r) { return r.pipv == 1 && r.pimv == 1; }
inline float DzPi(const R& r) { return std::max(std::fabs(r.pipz - r.vz), std::fabs(r.pimz - r.vz)); }

struct Win { double n[3] = {0, 0, 0}, var0 = 0; };
struct Sp {
    double h[3][kNb];
    double w[3][kNw];
    double cbW[kNw];
    Win win[2];       // 0 peak, 1 bump
};

// piecewise-linear CB density between the 10 MeV bin centres (constant outside)
// d(cb)/d(m) at mass m expressed as weights of the wide bins: m -> (j0, j1, u)
inline void Hat(double m, int& j0, int& j1, double& u) {
    const double t = (m - (kLo + 0.005)) / 0.010;
    if (t <= 0) { j0 = j1 = 0; u = 0; return; }
    if (t >= kNw - 1) { j0 = j1 = kNw - 1; u = 0; return; }
    j0 = (int)t; j1 = j0 + 1; u = t - j0;
}
// coefficients c_j of the CB integral over [a, b): CB = sum_j c_j cbW_j (exact, numerical)
std::vector<double> Coef(double a, double b, int ns = 200000) {
    std::vector<double> c(kNw, 0.);
    const double h = (b - a) / ns;
    for (int s = 0; s < ns; ++s) {
        const double m = a + (s + 0.5) * h;
        int j0, j1; double u; Hat(m, j0, j1, u);
        c[j0] += (1 - u) * h / 0.010; c[j1] += u * h / 0.010;
    }
    return c;
}

std::vector<Sp> Fill(const std::vector<R>* d, const std::vector<Sel>& sels) {
    const int ns = sels.size();
    std::vector<Sp> out(ns);
    for (auto& s : out) { std::memset(&s, 0, sizeof(Sp)); }
    for (int k = 0; k < 3; ++k) {
        for (const R& r : d[k]) {
            for (int s = 0; s < ns; ++s) {
                if (!sels[s](r)) continue;
                const double m = r.m4;
                if (!(m == m)) continue;
                if (m >= kLo && m < kHi) { const int ib = (int)((m - kLo) / kBw); out[s].h[k][ib] += 1; }
                if (m >= kPkLo && m < kPkHi) out[s].win[0].n[k] += 1;
                if (m >= kBmLo && m < kBmHi) out[s].win[1].n[k] += 1;
            }
        }
        // event-group variance of the unlike-sign window counts
        if (k == 0) {
            for (int wdx = 0; wdx < 2; ++wdx) {
                const double a = wdx == 0 ? kPkLo : kBmLo, b = wdx == 0 ? kPkHi : kBmHi;
                for (int s = 0; s < ns; ++s) {
                    int g = -1; double c = 0, var = 0;
                    for (const R& r : d[0]) {
                        if (r.grp != g) { var += c * c; c = 0; g = r.grp; }
                        if (sels[s](r) && r.m4 >= a && r.m4 < b) c += 1;
                    }
                    var += c * c;
                    out[s].win[wdx].var0 = var;
                }
            }
        }
    }
    for (auto& s : out) {
        for (int k = 0; k < 3; ++k) for (int i = 0; i < kNb; ++i) s.w[k][i / 4] += s.h[k][i];
        for (int j = 0; j < kNw; ++j) s.cbW[j] = 2.0 * std::sqrt(s.w[1][j] * s.w[2][j]);
    }
    return out;
}

void CBwin(const Sp& s, double a, double b, double& cb, double& ecb) {
    const std::vector<double> c = Coef(a, b);
    cb = 0; double v = 0;
    for (int j = 0; j < kNw; ++j) { cb += c[j] * s.cbW[j]; v += c[j] * c[j] * (s.w[1][j] + s.w[2][j]); }
    ecb = std::sqrt(v);
}

// CB per fine bin (exact integral over the bin) and its diagonal variance, scaled by f
void CBfine(const Sp& s, double f, std::vector<double>& cb, std::vector<double>& var) {
    static std::vector<std::vector<double>> cf;
    if (cf.empty()) for (int i = 0; i < kNb; ++i) cf.push_back(Coef(kLo + i * kBw, kLo + (i + 1) * kBw, 2000));
    cb.assign(kNb, 0.); var.assign(kNb, 0.);
    for (int i = 0; i < kNb; ++i) {
        const std::vector<double>& c = cf[i];
        for (int j = 0; j < kNw; ++j) if (c[j] != 0) { cb[i] += f * c[j] * s.cbW[j]; var[i] += f * f * c[j] * c[j] * (s.w[1][j] + s.w[2][j]); }
    }
}

const char* kFormula = "[0]*0.0025*TMath::Gaus(x,0.5439,0.0054,1) + [1] + [2]*((x-0.56625)/0.04375) + [3]*((x-0.56625)/0.04375)^2";

// P_main-like chi2 fit (ROOT) on N - f CB, bins 0.5225-0.6100
void FitMain(const Sp& s, double f, double& S, double& eS, double& chi2, int& ndf) {
    std::vector<double> cb, var; CBfine(s, f, cb, var);
    const int i0 = (int)std::lround((0.5225 - kLo) / kBw), i1 = (int)std::lround((0.6100 - kLo) / kBw);
    TH1D h("hmain", "", i1 - i0, 0.5225, 0.6100); h.SetDirectory(nullptr);
    for (int i = i0; i < i1; ++i) { h.SetBinContent(i - i0 + 1, s.h[0][i] - cb[i]); h.SetBinError(i - i0 + 1, std::sqrt(std::max(1.0, s.h[0][i] + var[i]))); }
    TF1 fn("fmain", kFormula, 0.5225, 0.6100);
    fn.SetParameters(1000, h.GetBinContent(1), 0, 0);
    TFitResultPtr r = h.Fit(&fn, "Q N 0 I S");
    S = fn.GetParameter(0); eS = fn.GetParError(0); chi2 = r->Chi2(); ndf = r->Ndf();
}
// Poisson likelihood fit (ROOT) of a like-sign spectrum (k = 1, 2, 3 = sum)
void FitLS(const Sp& s, int k, double& S, double& eS) {
    const int i0 = (int)std::lround((0.5225 - kLo) / kBw), i1 = (int)std::lround((0.6100 - kLo) / kBw);
    TH1D h("hls", "", i1 - i0, 0.5225, 0.6100); h.SetDirectory(nullptr);
    for (int i = i0; i < i1; ++i) { const double y = k == 3 ? s.h[1][i] + s.h[2][i] : s.h[k][i]; h.SetBinContent(i - i0 + 1, y); h.SetBinError(i - i0 + 1, std::sqrt(std::max(1.0, y))); }
    TF1 fn("fls", kFormula, 0.5225, 0.6100);
    fn.SetParameters(0, h.GetBinContent(18), 0, 0);
    h.Fit(&fn, "Q N 0 L I S");
    S = fn.GetParameter(0); eS = fn.GetParError(0);
}

} // namespace vcb

void verify_cb_counts() {
    using namespace vcb;
    std::vector<R> d[3];
    d[0] = Read("epem"); d[1] = Read("epep"); d[2] = Read("emem");

    std::vector<std::string> nm;
    std::vector<Sel> sv;
    auto add = [&](const char* n, Sel s) { nm.push_back(n); sv.push_back(s); };
    add("FAT only", [](const R& r) { return Fat(r); });
    add("pre", [](const R& r) { return Pre(r); });
    add("pre + OA>2", [](const R& r) { return Pre(r) && r.oa > 2.f; });
    add("pre + vertex", [](const R& r) { return Pre(r) && Vtx(r); });
    add("DIRECT", [](const R& r) { return Pre(r) && DzPi(r) < 40.f; });
    add("PEAK", [](const R& r) { return Pre(r) && DzPi(r) < 40.f && r.mm4 < 2.85f; });
    add("peak nominal", [](const R& r) { return Pre(r) && r.oalab < 55.f && r.mpipi > 0.30f; });
    add("peak nominal + vertex", [](const R& r) { return Pre(r) && r.oalab < 55.f && r.mpipi > 0.30f && Vtx(r); });
    add("bump nominal", [](const R& r) { return Pre(r) && r.oalab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f; });
    add("bump nominal + vertex", [](const R& r) { return Pre(r) && r.oalab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f && Vtx(r); });
    add("X17 (pre+Mee+vtx)", [](const R& r) { return Pre(r) && r.mee >= 0.0155f && r.mee < 0.0245f && Vtx(r); });
    add("pre, isBest==1", [](const R& r) { return Pre(r) && r.isBest == 1; });
    const std::vector<Sp> sp = Fill(d, sv);

    printf("\n1. WINDOW NUMBERS (own code). CB integrated exactly over the window; err(N-CB) = sqrt(event-group var + var CB)\n");
    printf("%-24s | %-4s | %8s %6s %6s %5s | %8s %6s | %6s | %6s | %9s %5s\n", "selection", "win", "N", "N++", "N--", "r", "CB", "eCB", "CB/N", "AM/GM", "N-CB", "err");
    for (size_t s = 0; s < sp.size(); ++s)
        for (int w = 0; w < 2; ++w) {
            const Win& x = sp[s].win[w];
            double cb, ecb; CBwin(sp[s], w == 0 ? kPkLo : kBmLo, w == 0 ? kPkHi : kBmHi, cb, ecb);
            printf("%-24s | %-4s | %8.0f %6.0f %6.0f %5.2f | %8.1f %6.1f | %6.4f | %6.3f | %9.0f %5.0f\n", nm[s].c_str(), w == 0 ? "peak" : "bump",
                   x.n[0], x.n[1], x.n[2], x.n[2] / x.n[1], cb, ecb, cb / x.n[0], (x.n[1] + x.n[2]) / (2 * std::sqrt(x.n[1] * x.n[2])),
                   x.n[0] - cb, std::sqrt(x.var0 + ecb * ecb));
        }
    // kept fractions of the vertex cuts
    {
        double cbF, cbP, cbV, e;
        for (int w = 0; w < 2; ++w) {
            const double a = w == 0 ? kPkLo : kBmLo, b = w == 0 ? kPkHi : kBmHi;
            CBwin(sp[0], a, b, cbF, e); CBwin(sp[1], a, b, cbP, e); CBwin(sp[3], a, b, cbV, e);
            printf("kept fractions, %s window: vz>-160 on FAT: N %.3f CB %.3f ; vertex on pre: N %.3f CB %.3f\n", w == 0 ? "peak" : "bump",
                   sp[1].win[w].n[0] / sp[0].win[w].n[0], cbP / cbF, sp[3].win[w].n[0] / sp[1].win[w].n[0], cbV / cbP);
        }
        double cbD; CBwin(sp[4], kPkLo, kPkHi, cbD, e); CBwin(sp[1], kPkLo, kPkHi, cbP, e);
        printf("dz_pi < 40 on pre, peak window: N %.3f CB %.3f\n", sp[4].win[0].n[0] / sp[1].win[0].n[0], cbD / cbP);
    }
    // whole sample like-sign ratio
    {
        double n[3] = {0, 0, 0};
        for (int k = 0; k < 3; ++k) for (const R& r : d[k]) if (Pre(r)) n[k] += 1;
        printf("pre, any M4: rows us %.0f pp %.0f mm %.0f, N--/N++ = %.3f\n", n[0], n[1], n[2], n[2] / n[1]);
    }

    printf("\n2. NULL TEST: ROOT Poisson likelihood fit (L I) of the frozen eta shape + pol2 on like-sign spectra, 0.5225-0.6100\n");
    printf("3. P_main-like ROOT chi2 fit (I) of N - CB, errors sqrt(N + var CB) without CB covariance; S for CB x 0.97, 1.00, 1.03\n");
    for (size_t s = 0; s < sp.size(); ++s) {
        double S, eS, c2, Sp_, eSp, Sm, eSm, Ss, eSs, S97, S103, e1; int ndf;
        FitMain(sp[s], 1.0, S, eS, c2, ndf);
        FitMain(sp[s], 0.97, S97, e1, c2, ndf);
        FitMain(sp[s], 1.03, S103, e1, c2, ndf);
        FitMain(sp[s], 1.0, S, eS, c2, ndf);
        FitLS(sp[s], 1, Sp_, eSp); FitLS(sp[s], 2, Sm, eSm); FitLS(sp[s], 3, Ss, eSs);
        printf("%-24s | S(N-CB) %7.1f +- %5.1f (chi2 %.1f/%d) | dS/S x0.97 %+.4f x1.03 %+.4f | S(++) %6.1f +- %5.1f | S(--) %6.1f +- %5.1f | S(sum) %6.1f +- %5.1f\n",
               nm[s].c_str(), S, eS, c2, ndf, (S97 - S) / S, (S103 - S) / S, Sp_, eSp, Sm, eSm, Ss, eSs);
    }

    // ---- 4. class-3 isBest test ----
    printf("\n4. LEPTON CLASSES (join by file, run, event, vertex x/y/z), selection pre, bump window\n");
    typedef std::tuple<float, float, float, uint32_t, uint32_t, uint32_t> Key;
    auto bits = [](float f) { uint32_t u; std::memcpy(&u, &f, 4); return u; };
    auto key = [&](const R& r) { return Key(r.file, r.run, r.ev, bits(r.vx), bits(r.vy), bits(r.vz)); };
    std::map<Key, std::pair<int, int>> mp;   // e+e- event -> (n e+, n e-)
    {
        size_t i = 0;
        while (i < d[0].size()) {
            size_t j = i;
            std::set<std::pair<float, float>> p, m;
            while (j < d[0].size() && d[0][j].grp == d[0][i].grp) { p.insert({d[0][j].epth, d[0][j].epph}); m.insert({d[0][j].emth, d[0][j].emph}); ++j; }
            mp[key(d[0][i])] = std::make_pair((int)p.size(), (int)m.size());
            i = j;
        }
    }
    double N[2][3][3];   // [best][tree][class 2,3,4]
    std::memset(N, 0, sizeof(N));
    long nmatch[3] = {0, 0, 0}, ngrp[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k) {
        size_t i = 0;
        while (i < d[k].size()) {
            size_t j = i;
            std::set<std::pair<float, float>> l1, l2, all;
            while (j < d[k].size() && d[k][j].grp == d[k][i].grp) {
                all.insert({d[k][j].epth, d[k][j].epph}); all.insert({d[k][j].emth, d[k][j].emph});
                l1.insert({d[k][j].epth, d[k][j].epph}); l2.insert({d[k][j].emth, d[k][j].emph});
                ++j;
            }
            int cls;
            ngrp[k]++;
            if (k == 0) cls = (int)(l1.size() + l2.size());
            else {
                int other = 0;
                auto it = mp.find(key(d[k][i]));
                if (it != mp.end()) { nmatch[k]++; other = (k == 1) ? it->second.second : it->second.first; }
                cls = (int)all.size() + other;
            }
            const int c = std::min(4, cls) - 2;
            for (size_t q = i; q < j; ++q) {
                const R& r = d[k][q];
                if (!Pre(r) || !(r.m4 >= kBmLo && r.m4 < kBmHi) || c < 0) continue;
                N[0][k][c] += 1;
                if (r.isBest == 1) N[1][k][c] += 1;
            }
            i = j;
        }
    }
    printf("   e+e- events %zu; e+e+ groups %ld matched %ld; e-e- groups %ld matched %ld\n", mp.size(), ngrp[1], nmatch[1], ngrp[2], nmatch[2]);
    double gs = 0, ge = 0, gsb = 0;
    const double eps[3] = {1.0, 0.5, 0.25};
    for (int c = 0; c < 3; ++c) {
        const double ga = 2 * std::sqrt(N[0][1][c] * N[0][2][c]), gb = 2 * std::sqrt(N[1][1][c] * N[1][2][c]);
        gs += ga; ge += eps[c] * ga; gsb += gb;
        printf("   class %d: N %8.0f N++ %6.0f N-- %6.0f GM %8.1f | isBest==1: N %8.0f GM %8.1f | eps_us %.3f eps_ls %.3f ratio %.3f\n", c + 2,
               N[0][0][c], N[0][1][c], N[0][2][c], ga, N[1][0][c], gb, N[1][0][c] / N[0][0][c], gb / ga, (N[1][0][c] / N[0][0][c]) / (gb / ga));
    }
    printf("   class shares of GM: %.3f %.3f %.3f ; X = sum GM / sum eps GM = %.3f\n", 2 * std::sqrt(N[0][1][0] * N[0][2][0]) / gs,
           2 * std::sqrt(N[0][1][1] * N[0][2][1]) / gs, 2 * std::sqrt(N[0][1][2] * N[0][2][2]) / gs, gs / ge);
}
