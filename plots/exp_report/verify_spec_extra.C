// ============================================================================
// verify_spec_extra.C  --  independent checks of the kinematic numbers of the
//                          chapter "spectra and kinematics" (prefix spec)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Written by the verifier of the spec chapter; does not use ScanEngine.h or
// spec_common.h; reads the FAT outputs (scan_nt) directly.
//   1. selection pre (cor_ok == 1, isBest >= 0, vz > -160 mm) and pre + MM4 <
//      2.85 GeV/c2: N - CB in the data peak / bump windows (CB = 2 sqrt(N++ N--)
//      per 10 MeV M4 bin, piecewise linear density between bin centres,
//      integrated over the window);
//   2. distribution of the window CB over another variable x (own code of the
//      "M4-binned" method: a like-sign row of 10 MeV bin j gets the weight
//      c_j sqrt(N_other_j / N_own_j), c_j = integral of the linear hat of bin j
//      over the window): fractions of N - CB with M(ee) < 10, < 20, > 50 MeV/c2,
//      OA < 2, < 6 deg, and means of N - CB for M(ee), OA, M(pipi), oa_lab,
//      MM4, pT, y_cm, eta-rest angle, dz_pi (all M4, peak window, bump window);
//      also the "subset" alternative: CB recomputed from the like-sign rows of
//      the subset only (2 sqrt(N++ N--) of the subset, per 10 MeV bin);
//   3. vz of cor_ok, isBest >= 0 rows: most populated 5 mm bin, fraction in
//      -200 <= vz < -160 mm, per sample.
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/output_pippim{epem,epep,emem}_exp.root (scan_nt)
// Outputs: printout only
// Run (from /home/przygoda/HADES/FAT/FAT; about 30 s):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_spec_verify", kTRUE)' 'plots/exp_report/verify_spec_extra.C+'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

namespace vse {
struct R { float best, cor, vz, m4, mee, oa, mpipi, oalab, mm4, pt4, y4, oaeta, dzpi; };
const char* const kDir = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530";
const char* const kFile[3] = {"output_pippimepem_exp.root", "output_pippimepep_exp.root", "output_pippimemem_exp.root"};

std::vector<R> Read(int k) {
    std::vector<R> v;
    TFile* f = TFile::Open(Form("%s/%s", kDir, kFile[k]));
    TTree* t = (TTree*)f->Get("scan_nt");
    const char* nm[] = {"isBest", "cor_ok", "eVertReco_z", "m_pippimepem_cor", "m_ee_cor", "oa", "m_pippim_cor",
                        "oa_pippim_epem_lab_cor", "mm_pippimepem_cor", "pt4_cor", "y4_cms_cor", "oa_pippim_epem_eta_rest_cor", "pip_z", "pim_z"};
    const int n = 14;
    Float_t x[14];
    t->SetBranchStatus("*", 0);
    for (int i = 0; i < n; ++i) { t->SetBranchStatus(nm[i], 1); t->SetBranchAddress(nm[i], &x[i]); }
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        const float a = std::fabs(x[12] - x[2]), b = std::fabs(x[13] - x[2]);
        v.push_back({x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], a > b ? a : b});
    }
    f->Close();
    return v;
}
// coefficients c_j of the 10 MeV bins (0.25 + 0.01 j) for the window [a, b): integral of the linear hat
std::vector<double> Coef(double a, double b) {
    const int nb = 300; const double lo = 0.25, w = 0.01;
    std::vector<double> c(nb, 0.);
    const int nst = (int)std::ceil((b - a) / 1e-6); const double h = (b - a) / nst;
    for (int q = 0; q < nst; ++q) {
        const double xx = a + (q + 0.5) * h;
        const double t = (xx - lo) / w - 0.5;
        int j0, j1; double f;
        if (t <= 0) { j0 = j1 = 0; f = 0; } else if (t >= nb - 1) { j0 = j1 = nb - 1; f = 0; } else { j0 = (int)t; j1 = j0 + 1; f = t - j0; }
        c[j0] += (1 - f) * h / w; c[j1] += f * h / w;
    }
    return c;
}
}  // namespace vse

void verify_spec_extra() {
    using namespace vse;
    std::vector<R> d[3];
    for (int k = 0; k < 3; ++k) d[k] = Read(k);
    auto pre = [](const R& r) { return r.cor == 1.f && r.best > -0.5f && r.vz > -160.f; };
    const double pkLo = 0.5439 - 1.9 * 0.0054, pkHi = 0.5439 + 1.9 * 0.0054;
    const double wl[3][2] = {{0.25, 3.25}, {pkLo, pkHi}, {0.455, 0.5225}};
    const char* wn[3] = {"all", "peak", "bump"};
    std::vector<double> cj[3];
    for (int w = 0; w < 3; ++w) cj[w] = Coef(wl[w][0], wl[w][1]);

    // generic: N - CB in window w for selection s, with CB from the like-sign rows of s (subset method)
    auto nmcb = [&](int w, const std::function<bool(const R&)>& s, double* cbOut) {
        std::vector<double> np(300, 0.), nm(300, 0.); double n = 0;
        for (int k = 0; k < 3; ++k)
            for (const R& r : d[k]) {
                if (!s(r) || !(r.m4 == r.m4)) continue;
                if (k == 0) { if (r.m4 >= wl[w][0] && r.m4 < wl[w][1]) n += 1; continue; }
                const int j = (int)((r.m4 - 0.25) / 0.01);
                if (j < 0 || j >= 300) continue;
                (k == 1 ? np : nm)[j] += 1;
            }
        double cb = 0;
        for (int j = 0; j < 300; ++j) cb += cj[w][j] * 2 * std::sqrt(np[j] * nm[j]);
        if (cbOut) *cbOut = cb;
        return n - cb;
    };
    // 1. MM4 < 2.85 on pre
    for (int w = 1; w < 3; ++w) {
        const double a = nmcb(w, pre, nullptr);
        const double b = nmcb(w, [&](const R& r) { return pre(r) && r.mm4 < 2.85f; }, nullptr);
        printf("window %s: N-CB pre %.0f, pre+MM4<2.85 %.0f, removed fraction %.4f\n", wn[w], a, b, 1 - b / a);
    }
    // 2. distributed CB (M4-binned) per window
    std::vector<double> np(300, 0.), nm(300, 0.);
    for (int k = 1; k < 3; ++k)
        for (const R& r : d[k]) {
            if (!pre(r)) continue;
            const int j = (int)((r.m4 - 0.25) / 0.01);
            if (j >= 0 && j < 300) (k == 1 ? np : nm)[j] += 1;
        }
    typedef std::function<double(const R&)> V;
    struct Q { const char* name; V x; };
    std::vector<Q> qs = {{"mee[MeV]", [](const R& r) { return 1000. * r.mee; }}, {"oa", [](const R& r) { return (double)r.oa; }},
                         {"mpipi", [](const R& r) { return (double)r.mpipi; }}, {"oalab", [](const R& r) { return (double)r.oalab; }},
                         {"mm4", [](const R& r) { return (double)r.mm4; }}, {"pt4[GeV]", [](const R& r) { return 1e-3 * r.pt4; }},
                         {"y4", [](const R& r) { return (double)r.y4; }},
                         {"oaeta", [](const R& r) { return (double)r.oaeta; }}, {"dzpi", [](const R& r) { return (double)r.dzpi; }}};
    for (int w = 0; w < 3; ++w) {
        // sums: N, CB, and N*x, CB*x, plus fractions
        double sN = 0, sC = 0;
        std::vector<double> sNx(qs.size(), 0.), sCx(qs.size(), 0.);
        double fN[5] = {0, 0, 0, 0, 0}, fC[5] = {0, 0, 0, 0, 0};
        auto cls = [](const R& r, int i) {
            const double m = 1000. * r.mee;
            switch (i) { case 0: return m < 10; case 1: return m < 20; case 2: return m > 50; case 3: return r.oa < 2.f; default: return r.oa < 6.f; }
        };
        for (int k = 0; k < 3; ++k)
            for (const R& r : d[k]) {
                if (!pre(r)) continue;
                double wt;
                if (k == 0) { if (!(r.m4 >= wl[w][0] && r.m4 < wl[w][1])) continue; wt = 1; }
                else {
                    const int j = (int)((r.m4 - 0.25) / 0.01);
                    if (j < 0 || j >= 300 || cj[w][j] == 0 || np[j] == 0 || nm[j] == 0) continue;
                    wt = cj[w][j] * (k == 1 ? std::sqrt(nm[j] / np[j]) : std::sqrt(np[j] / nm[j]));
                }
                double& S = (k == 0) ? sN : sC;
                S += wt;
                for (size_t q = 0; q < qs.size(); ++q) (k == 0 ? sNx : sCx)[q] += wt * qs[q].x(r);
                for (int i = 0; i < 5; ++i) if (cls(r, i)) (k == 0 ? fN : fC)[i] += wt;
            }
        printf("\nwindow %s: N %.0f CB(distributed) %.1f\n", wn[w], sN, sC);
        printf("  fractions of N - CB: mee<10 %.3f  mee<20 %.3f  mee>50 %.3f  OA<2 %.3f  OA<6 %.3f\n", (fN[0] - fC[0]) / (sN - sC),
               (fN[1] - fC[1]) / (sN - sC), (fN[2] - fC[2]) / (sN - sC), (fN[3] - fC[3]) / (sN - sC), (fN[4] - fC[4]) / (sN - sC));
        printf("  means of N - CB:");
        for (size_t q = 0; q < qs.size(); ++q) printf(" %s %.4g", qs[q].name, (sNx[q] - sCx[q]) / (sN - sC));
        printf("\n  means of N     :");
        for (size_t q = 0; q < qs.size(); ++q) printf(" %s %.4g", qs[q].name, sNx[q] / sN);
        printf("\n  means of CB    :");
        for (size_t q = 0; q < qs.size(); ++q) printf(" %s %.4g", qs[q].name, sCx[q] / sC);
        printf("\n");
        // subset alternative for the class fractions
        if (w > 0) {
            const double tot = nmcb(w, pre, nullptr);
            printf("  subset-CB fractions:");
            for (int i = 0; i < 5; ++i) printf(" %.3f", nmcb(w, [&](const R& r) { return pre(r) && cls(r, i); }, nullptr) / tot);
            printf("   (mee<10, mee<20, mee>50, OA<2, OA<6)\n");
        }
    }
    // 3. vz shapes
    for (int k = 0; k < 3; ++k) {
        TH1D h("h", "", 120, -500, 100); h.SetDirectory(nullptr);
        double n = 0, a = 0;
        for (const R& r : d[k]) {
            if (!(r.cor == 1.f && r.best > -0.5f)) continue;
            n += 1; h.Fill(r.vz);
            if (r.vz >= -200.f && r.vz < -160.f) a += 1;
        }
        const int mb = h.GetMaximumBin();
        printf("sample %d: vz max bin %.0f..%.0f mm, fraction in [-200,-160) %.4f, rows above 100 mm %.0f (of %.0f)\n", k,
               h.GetXaxis()->GetBinLowEdge(mb), h.GetXaxis()->GetBinUpEdge(mb), a / n, h.GetBinContent(121), n);
    }
}
