// ============================================================================
// verify_spec_counts.C  --  independent re-derivation of the key numbers of the
//                           "spectra and kinematics" chapter (prefix spec)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Written by the verifier of the spec chapter. It does NOT use ScanEngine.h,
// spec_common.h or the compact file: it reads the three FAT outputs (tree
// scan_nt) directly and recomputes, with its own code,
//   - row totals of pre (cor_ok == 1 && isBest >= 0 && vz > -160 mm);
//   - N, N++, N-- in the data peak window 0.53364-0.55416 and the data bump
//     window 0.4550-0.5225 GeV/c2 for the 8 selections x {isBest >= 0,
//     isBest == 1} x {no OA cut, OA > 2 deg}; event-group error of N
//     (group = adjacent rows with equal file_id, runnumber, event);
//   - CB = 2 sqrt(N++ N--) per 10 MeV M4 bin, as a piecewise linear density
//     between the 10 MeV bin centres (flat beyond the outer centres), integrated
//     over the window (fine numerical integration; not the fine-bin sum of
//     ScanEngine); Poisson variance of the CB (no event-group term);
//   - visible yields (N - CB) x 64/63 / 6079.2 nb^-1;
//   - an own P_main: TH1 of N - CB in 2.5 MeV bins (0.25-3.25), errors
//     sqrt(N + var CB), ROOT/Minuit chi2 fit of S x (Gaussian integral over the
//     bin, 0.5439 GeV/c2, 5.4 MeV) + pol2 in bins with centres in 0.5225-0.6100;
//     the same with CB x 1.03 and x 0.97;
//   - row pass fractions of the vertex cuts (vz > -160 on cor_ok, both pions in
//     the vertex fit on pre, dz_pi < 40 mm on pre) in the two windows.
// Selections as defined for the report:
//   pre = cor_ok == 1 && isBest >= 0 && vz > -160 (FAT has vz > -500 and PT3)
//   +vtx = pip_isUsedVertex == 1 && pim_isUsedVertex == 1
//   peak nominal = pre + oa_lab < 55 + M(pipi) > 0.30; bump nominal = pre +
//   oa_lab < 40 + 0.30 < M(pipi) < 0.41; PEAK = pre + dz_pi < 40 + MM4 < 2.85;
//   DIRECT = pre + dz_pi < 40; dz_pi = max(|pip_z - vz|, |pim_z - vz|).
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/output_pippim{epem,epep,emem}_exp.root (scan_nt)
// Outputs: printout only (the verifier copies the numbers into spec_results.md)
// Run (from /home/przygoda/HADES/FAT/FAT; about 1 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_spec_verify", kTRUE)' 'plots/exp_report/verify_spec_counts.C+'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "Math/MinimizerOptions.h"

namespace vsc {

struct R {
    int grp;
    float best, cor, vz, m4, mpipi, mee, oa, oalab, mm4, dzpi, pipv, pimv;
};

const char* const kDir = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530";
const char* const kFile[3] = {"output_pippimepem_exp.root", "output_pippimepep_exp.root", "output_pippimemem_exp.root"};
const double kPkLo = 0.5439 - 1.9 * 0.0054, kPkHi = 0.5439 + 1.9 * 0.0054;   // 0.53364 - 0.55416
const double kBpLo = 0.4550, kBpHi = 0.5225;
const double kL = 6079.2, kPT3 = 64.0 / 63.0;

std::vector<R> ReadSample(int k, std::map<float, int>& vtxVals) {
    std::vector<R> v;
    TFile* f = TFile::Open(Form("%s/%s", kDir, kFile[k]));
    TTree* t = (TTree*)f->Get("scan_nt");
    const char* nm[] = {"file_id", "runnumber", "event", "isBest", "cor_ok", "eVertReco_z", "m_pippimepem_cor", "m_pippim_cor",
                        "m_ee_cor", "oa", "oa_pippim_epem_lab_cor", "mm_pippimepem_cor", "pip_z", "pim_z", "pip_isUsedVertex",
                        "pim_isUsedVertex", "trigbit"};
    const int n = sizeof(nm) / sizeof(nm[0]);
    std::vector<Float_t> x(n);
    t->SetBranchStatus("*", 0);
    for (int i = 0; i < n; ++i) { t->SetBranchStatus(nm[i], 1); t->SetBranchAddress(nm[i], &x[i]); }
    const Long64_t ne = t->GetEntries();
    v.reserve(ne);
    float pf = -1, pr = -1, pe = -1; int g = -1;
    long nBadTrig = 0, nVzBad = 0;
    for (Long64_t i = 0; i < ne; ++i) {
        t->GetEntry(i);
        if (x[0] != pf || x[1] != pr || x[2] != pe) { ++g; pf = x[0]; pr = x[1]; pe = x[2]; }
        if (x[16] != 8192.f) ++nBadTrig;
        if (!(x[5] > -500.f)) ++nVzBad;
        vtxVals[x[14]]++;
        R r;
        r.grp = g; r.best = x[3]; r.cor = x[4]; r.vz = x[5]; r.m4 = x[6]; r.mpipi = x[7]; r.mee = x[8]; r.oa = x[9];
        r.oalab = x[10]; r.mm4 = x[11];
        const float a = std::fabs(x[12] - x[5]), b = std::fabs(x[13] - x[5]);
        r.dzpi = a > b ? a : b;
        r.pipv = x[14]; r.pimv = x[15];
        v.push_back(r);
    }
    printf("sample %d: %lld rows, %d groups, trigbit != 8192: %ld, vz <= -500: %ld\n", k, ne, g + 1, nBadTrig, nVzBad);
    f->Close();
    return v;
}

typedef std::function<bool(const R&)> Sel;

// CB integral over [a, b) from the 10 MeV like-sign counts (piecewise linear density)
double CBint(const TH1D* pp10, const TH1D* mm10, double a, double b, double* var) {
    const int nb = pp10->GetNbinsX();
    std::vector<double> cb(nb), cv(nb), c(nb, 0.);
    for (int j = 0; j < nb; ++j) {
        const double p = pp10->GetBinContent(j + 1), m = mm10->GetBinContent(j + 1);
        cb[j] = 2 * std::sqrt(p * m);
    }
    const double lo = pp10->GetXaxis()->GetXmin(), w = pp10->GetBinWidth(1);
    // density at x: linear between centres; coefficient bookkeeping for the variance
    const int nst = (int)std::ceil((b - a) / 1e-6);
    const double h = (b - a) / nst;
    double s = 0;
    for (int q = 0; q < nst; ++q) {
        const double xx = a + (q + 0.5) * h;
        double t = (xx - lo) / w - 0.5;
        int j0, j1; double f;
        if (t <= 0) { j0 = j1 = 0; f = 0; }
        else if (t >= nb - 1) { j0 = j1 = nb - 1; f = 0; }
        else { j0 = (int)t; j1 = j0 + 1; f = t - j0; }
        c[j0] += (1 - f) * h / w; c[j1] += f * h / w;
    }
    double v = 0;
    for (int j = 0; j < nb; ++j) {
        s += c[j] * cb[j];
        const double p = pp10->GetBinContent(j + 1), m = mm10->GetBinContent(j + 1);
        if (p > 0 && m > 0) v += c[j] * c[j] * (p + m);   // dCB/dp = sqrt(m/p): var = m/p * p + p/m * m
    }
    if (var) *var = v;
    return s;
}

struct Res { double n, nvar, npp, nmm, cb, cbvar; };

}  // namespace vsc

void verify_spec_counts() {
    using namespace vsc;
    TH1::AddDirectory(kFALSE);
    TStopwatch sw;
    std::map<float, int> vv[3];
    std::vector<R> d[3];
    for (int k = 0; k < 3; ++k) d[k] = ReadSample(k, vv[k]);
    for (int k = 0; k < 3; ++k) {
        printf("sample %d pip_isUsedVertex values:", k);
        for (auto& e : vv[k]) printf(" %g:%d", e.first, e.second);
        printf("\n");
    }
    auto pre = [](const R& r) { return r.cor == 1.f && r.best > -0.5f && r.vz > -160.f; };
    auto vtx = [](const R& r) { return r.pipv == 1.f && r.pimv == 1.f; };
    auto pk = [](const R& r) { return r.oalab < 55.f && r.mpipi > 0.30f; };
    auto bp = [](const R& r) { return r.oalab < 40.f && r.mpipi > 0.30f && r.mpipi < 0.41f; };
    std::vector<std::pair<std::string, Sel>> base = {
        {"pre", [=](const R& r) { return pre(r); }},
        {"pre+vtx", [=](const R& r) { return pre(r) && vtx(r); }},
        {"peaknom", [=](const R& r) { return pre(r) && pk(r); }},
        {"peaknom+vtx", [=](const R& r) { return pre(r) && pk(r) && vtx(r); }},
        {"bumpnom", [=](const R& r) { return pre(r) && bp(r); }},
        {"bumpnom+vtx", [=](const R& r) { return pre(r) && bp(r) && vtx(r); }},
        {"PEAK", [=](const R& r) { return pre(r) && r.dzpi < 40.f && r.mm4 < 2.85f; }},
        {"DIRECT", [=](const R& r) { return pre(r) && r.dzpi < 40.f; }},
    };
    const char* vn[4] = {"best>=0", "best>=0,OA>2", "best==1", "best==1,OA>2"};

    // totals of pre
    for (int k = 0; k < 3; ++k) {
        long n = 0; for (const R& r : d[k]) if (pre(r)) ++n;
        printf("pre total rows sample %d: %ld\n", k, n);
    }

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    printf("\n%-12s %-13s | %9s %7s | %8s %6s %6s %8s | %16s %7s | %16s %8s | %7s %7s | %20s %6s %6s\n", "selection", "variant", "N_pk",
           "eN", "CB_pk", "N++", "N--", "CB/N", "N-CB peak", "Y_pk", "N-CB bump", "Y_bp", "CB/N bp", "", "S (own fit)", "dS+3%", "dS-3%");
    for (size_t b = 0; b < base.size(); ++b) {
        for (int v = 0; v < 4; ++v) {
            const Sel s0 = base[b].second;
            Sel s = [=](const R& r) {
                if (!s0(r)) return false;
                if ((v == 1 || v == 3) && !(r.oa > 2.f)) return false;
                if (v >= 2 && !(r.best == 1.f)) return false;
                return true;
            };
            TH1D h[3] = {TH1D("h0", "", 1200, 0.25, 3.25), TH1D("h1", "", 1200, 0.25, 3.25), TH1D("h2", "", 1200, 0.25, 3.25)};
            Res win[2];
            for (int w = 0; w < 2; ++w) win[w] = {0, 0, 0, 0, 0, 0};
            for (int k = 0; k < 3; ++k) {
                h[k].SetDirectory(nullptr);
                int cur = -1; int kw[2] = {0, 0};
                for (const R& r : d[k]) {
                    if (r.grp != cur) {
                        if (k == 0) for (int w = 0; w < 2; ++w) { win[w].nvar += (double)kw[w] * kw[w]; kw[w] = 0; }
                        cur = r.grp;
                    }
                    if (!s(r)) continue;
                    if (!(r.m4 == r.m4)) continue;
                    h[k].Fill(r.m4);
                    const double lo[2] = {kPkLo, kBpLo}, hi[2] = {kPkHi, kBpHi};
                    for (int w = 0; w < 2; ++w)
                        if (r.m4 >= lo[w] && r.m4 < hi[w]) {
                            if (k == 0) { win[w].n += 1; kw[w]++; }
                            else if (k == 1) win[w].npp += 1;
                            else win[w].nmm += 1;
                        }
                }
                if (k == 0) for (int w = 0; w < 2; ++w) { win[w].nvar += (double)kw[w] * kw[w]; kw[w] = 0; }
            }
            TH1D* pp10 = (TH1D*)h[1].Rebin(4, "pp10"); TH1D* mm10 = (TH1D*)h[2].Rebin(4, "mm10");
            pp10->SetDirectory(nullptr); mm10->SetDirectory(nullptr);
            win[0].cb = CBint(pp10, mm10, kPkLo, kPkHi, &win[0].cbvar);
            win[1].cb = CBint(pp10, mm10, kBpLo, kBpHi, &win[1].cbvar);
            // own P_main: N - CB per 2.5 MeV bin, CB integral over the bin
            double S[3] = {0, 0, 0}, eS = 0, chi2 = 0; int ndf = 0;
            for (int q = 0; q < 3; ++q) {
                const double fcb = q == 0 ? 1.0 : (q == 1 ? 1.03 : 0.97);
                TH1D hs("hs", "", 1200, 0.25, 3.25); hs.SetDirectory(nullptr);
                for (int i = 1; i <= 1200; ++i) {
                    const double x1 = hs.GetXaxis()->GetBinLowEdge(i), x2 = x1 + 0.0025;
                    if (x2 < 0.50 || x1 > 0.63) continue;
                    double cv = 0;
                    const double cb = fcb * CBint(pp10, mm10, x1, x2, &cv);
                    hs.SetBinContent(i, h[0].GetBinContent(i) - cb);
                    hs.SetBinError(i, std::sqrt(h[0].GetBinContent(i) + fcb * fcb * cv));
                }
                TF1 fn("fn", [](double* x, double* p) {
                    const double m = x[0], bw = 0.0025, mu = 0.5439, sg = 0.0054;
                    const double g = 0.5 * (TMath::Erf((m + 0.5 * bw - mu) / (std::sqrt(2.) * sg)) - TMath::Erf((m - 0.5 * bw - mu) / (std::sqrt(2.) * sg)));
                    const double t = (m - 0.56625) / 0.04375;
                    return p[0] * g + p[1] + p[2] * t + p[3] * t * t;
                }, 0.5225, 0.6100, 4);
                fn.SetParameters(4000, 300, -50, 0);
                hs.Fit(&fn, "QR0N");
                S[q] = fn.GetParameter(0);
                if (q == 0) { eS = fn.GetParError(0); chi2 = fn.GetChisquare(); ndf = fn.GetNDF(); }
            }
            const Res& p = win[0]; const Res& u = win[1];
            printf("%-12s %-13s | %9.0f %7.0f | %8.1f %6.0f %6.0f %8.4f | %9.0f +- %4.0f %7.3f | %9.0f +- %4.0f %8.3f | %7.4f %7s | %8.0f +- %4.0f (%4.1f/%2d) %+6.1f %+6.1f\n",
                   base[b].first.c_str(), vn[v], p.n, std::sqrt(p.nvar), p.cb, p.npp, p.nmm, p.cb / p.n, p.n - p.cb,
                   std::sqrt(p.nvar + p.cbvar), (p.n - p.cb) * kPT3 / kL, u.n - u.cb, std::sqrt(u.nvar + u.cbvar),
                   (u.n - u.cb) * kPT3 / kL, u.cb / u.n, "", S[0], eS, chi2, ndf, S[1] - S[0], S[2] - S[0]);
            delete pp10; delete mm10;
        }
    }

    // vertex-cut pass fractions in the windows
    printf("\nrow pass fractions (peak window / bump window)\n");
    const char* cn[3] = {"vz>-160 (base cor_ok, best>=0)", "pi pi in vertex fit (base pre)", "dz_pi<40 (base pre)"};
    for (int c = 0; c < 3; ++c) {
        printf("%-34s", cn[c]);
        for (int k = 0; k < 3; ++k) {
            double np[2] = {0, 0}, nt[2] = {0, 0};
            for (const R& r : d[k]) {
                bool bs, ps;
                if (c == 0) { bs = r.cor == 1.f && r.best > -0.5f; ps = r.vz > -160.f; }
                else { bs = pre(r); ps = (c == 1) ? vtx(r) : (r.dzpi < 40.f); }
                if (!bs) continue;
                const int w = (r.m4 >= kPkLo && r.m4 < kPkHi) ? 0 : ((r.m4 >= kBpLo && r.m4 < kBpHi) ? 1 : -1);
                if (w < 0) continue;
                nt[w] += 1; if (ps) np[w] += 1;
            }
            printf("  k%d %.4f / %.4f", k, np[0] / nt[0], np[1] / nt[1]);
        }
        printf("\n");
    }
    // fraction of vz <= -160 among cor_ok, isBest >= 0 rows
    for (int k = 0; k < 3; ++k) {
        double a = 0, t = 0;
        for (const R& r : d[k]) { if (!(r.cor == 1.f && r.best > -0.5f)) continue; t += 1; if (!(r.vz > -160.f)) a += 1; }
        printf("sample %d: fraction of cor_ok rows with vz <= -160: %.4f\n", k, a / t);
    }
    printf("done in %.0f s\n", sw.RealTime());
}
