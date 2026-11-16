// ============================================================================
// v0_verify.C  --  independent verification of WP0 (prefix v0)
// ============================================================================
// Recomputes histogram contents of hists_v1.root directly from the compact files,
// WITHOUT defs.h, hists_io.h or the WP0 copies of the simulation headers:
//   data : raw TTree reading of scan_v2/scan_compact.root (own selections, own
//          windows, own CB = 2 sqrt(N++ N--) in 10 MeV cells + linear interpolation +
//          f(cell | 50 MeV band)), event-group adjacency, dz_pi column, rows below 0.30
//   sim  : the ORIGINAL simulation-branch headers in FAT_simwt (read-only include),
//          streaming over std_us, lep_us, eta_us (no Load into memory), own mass
//          transfer, own MM4 at 4500 MeV, own event-group variance (global map)
// and compares every number with hists_v1.root (read-only).
// Output: docs/latex/tables_all/v0_verify.txt
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_v0", kTRUE)' 'plots/all_report/v0_verify.C+O'   (argument: 1 data, 2 sim, 4 file; default 7 = all)
// ============================================================================
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"

// ORIGINAL headers of the simulation study (read-only; they include their own SimEngine.h, ScanEngine.h, sim_truth.h)
#include "/home/przygoda/HADES/FAT/FAT_simwt/plots/eta_sim/sim_eta_physics.h"
#include "/home/przygoda/HADES/FAT/FAT_simwt/plots/eta_sim/sim_ghostvar.h"

namespace v0 {
const char* const kData = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root";
const char* const kSim = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/sim_compact.root";
const char* const kSimEta = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/sim_compact_eta.root";
const char* const kHists = "/hdd3/przygoda/HADES/FAT_out/all_report/hists_v1.root";
const char* const kOut = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/v0_verify.txt";
FILE* gOut = nullptr;
int gNfail = 0;

// window ids as in the "win" histogram (bin = id + 1)
enum { wB = 0, wP, wSBlo, wSBhi, wX17c, wCRU, wCRK, wCRL, wR0 = 8, wX17P = 20 };

void Cmp(const char* what, double mine, double file, double tol, const char* note = "") {
    const double d = std::fabs(mine - file);
    const bool ok = d <= tol;
    if (!ok) ++gNfail;
    char buf[1024];
    snprintf(buf, sizeof(buf), "%-62s mine %16.6f  file %16.6f  diff %10.3g  %s %s", what, mine, file, mine - file, ok ? "OK" : "MISMATCH", note);
    const std::string line(buf);
    printf("%s\n", line.c_str());
    fprintf(gOut, "%s\n", line.c_str());
}
void Note(const std::string& s) { printf("%s\n", s.c_str()); fprintf(gOut, "%s\n", s.c_str()); }

TH1* H(TFile* f, const char* p) {
    TH1* h = (TH1*)f->Get(p);
    if (!h) { Note(Form("MISSING in hists_v1.root: %s", p)); ++gNfail; }
    return h;
}

// ---------------------------------------------------------------------------- CB (own implementation)
const int NF = 320;            // fine bins 0.30-1.10, 2.5 MeV
const double FLO = 0.30, FW = 0.0025;
inline int Fine(double m) { if (!(m >= FLO && m < 1.10)) return -1; int i = (int)std::floor((m - FLO) / FW); return i > NF - 1 ? NF - 1 : i; }
// fine CB from fine like-sign counts: 10 MeV cells, 2 sqrt(pp mm), linear interpolation between cell centres
std::vector<double> FineCB(const std::vector<double>& pp, const std::vector<double>& mm, int reb = 4) {
    const int nw = NF / reb;
    std::vector<double> W(nw);
    for (int j = 0; j < nw; ++j) {
        double a = 0, b = 0;
        for (int q = 0; q < reb; ++q) { a += pp[j * reb + q]; b += mm[j * reb + q]; }
        W[j] = 2. * std::sqrt(a * b);
    }
    std::vector<double> cb(NF);
    const double cw = reb * FW;
    for (int i = 0; i < NF; ++i) {
        const double x = FLO + (i + 0.5) * FW;
        const double c0 = FLO + 0.5 * cw, cl = FLO + (nw - 0.5) * cw;
        double v;
        if (x <= c0) v = W[0];
        else if (x >= cl) v = W[nw - 1];
        else {
            const double t = (x - c0) / cw;
            const int j = (int)std::floor(t);
            const double f = t - j;
            v = (1 - f) * W[j] + f * W[j + 1];
        }
        cb[i] = v / reb;
    }
    return cb;
}
// CB of [a, b) with fractional fine bins, optional per-band factor fb[band]
double WinCB(const std::vector<double>& cb, double a, double b, const std::vector<double>* fb = nullptr) {
    double s = 0;
    for (int i = 0; i < NF; ++i) {
        const double x1 = FLO + i * FW, x2 = x1 + FW;
        const double o = std::min(x2, b) - std::max(x1, a);
        if (o <= 1e-12) continue;
        double fr = std::min(1., o / FW);
        if (fb) fr *= (*fb)[i / 20];
        s += fr * cb[i];
    }
    return s;
}

// ---------------------------------------------------------------------------- data
inline int SliceOf(double m) { if (!(m == m)) return -1; return m < 0.34 ? 0 : m < 0.37 ? 1 : m < 0.39 ? 2 : m < 0.43 ? 3 : 4; }

struct DSelAcc {
    // fine counts per sample (0 us, 1 pp, 2 mm)
    std::vector<double> fine[3];
    // band counts of like-sign rows: total with a valid slice, S5, X17 M(ee), F10 cell (slice, E) and total (any)
    std::vector<double> bandAll, bandSlice, bandS5, bandX17, bandCell[10];
    double win[3][21];
    double winVarGrp = 0;    // event-group variance of N in P (us)
    double f10[3];           // counts in the chosen F10 cell (S2E2, 0.540-0.5425)
    double x2d[3];           // X2D cell M4 0.5425-0.545, M(ee) 30-32 MeV
    double below030[3];
    double perIIA_P = 0, fold0_P = 0;
    DSelAcc() {
        for (int k = 0; k < 3; ++k) { fine[k].assign(NF, 0.); for (int w = 0; w < 21; ++w) win[k][w] = 0; f10[k] = x2d[k] = below030[k] = 0; }
        bandAll.assign(16, 0.); bandSlice.assign(16, 0.); bandS5.assign(16, 0.); bandX17.assign(16, 0.);
        for (int c = 0; c < 10; ++c) bandCell[c].assign(16, 0.);
    }
};

void DataPart(TFile* hf) {
    Note("=============================== DATA (raw TTree reading of scan_v2/scan_compact.root)");
    TFile* f = TFile::Open(kData, "READ");
    const char* trees[3] = {"us", "pp", "mm"};
    // selections: 0 S0, 1 pre, 2 none, 3 S0_corfb, 4 S0_best
    const char* stag[5] = {"S0", "pre", "none", "S0_corfb", "S0_best"};
    DSelAcc acc[5];
    long long nAdjViol[3] = {0, 0, 0}, nDzBad[3] = {0, 0, 0}, nDzChk[3] = {0, 0, 0}, nRows[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k) {
        TTree* t = (TTree*)f->Get(trees[k]);
        Int_t grp, fold, file_id, isBest, cor_ok, pipv, pimv;
        Float_t m4, mee, mpipi, vz, dz_pi, pip_z, pim_z, m4_fb, mee_fb, mpipi_fb;
        t->SetBranchStatus("*", 0);
        for (auto n : {"grp", "fold", "file_id", "isBest", "cor_ok", "pip_isUsedVertex", "pim_isUsedVertex", "m4", "mee", "mpipi", "vz",
                       "dz_pi", "pip_z", "pim_z", "m4_fb", "mee_fb", "mpipi_fb"})
            t->SetBranchStatus(n, 1);
        t->SetBranchAddress("grp", &grp); t->SetBranchAddress("fold", &fold); t->SetBranchAddress("file_id", &file_id);
        t->SetBranchAddress("isBest", &isBest); t->SetBranchAddress("cor_ok", &cor_ok);
        t->SetBranchAddress("pip_isUsedVertex", &pipv); t->SetBranchAddress("pim_isUsedVertex", &pimv);
        t->SetBranchAddress("m4", &m4); t->SetBranchAddress("mee", &mee); t->SetBranchAddress("mpipi", &mpipi); t->SetBranchAddress("vz", &vz);
        t->SetBranchAddress("dz_pi", &dz_pi); t->SetBranchAddress("pip_z", &pip_z); t->SetBranchAddress("pim_z", &pim_z);
        t->SetBranchAddress("m4_fb", &m4_fb); t->SetBranchAddress("mee_fb", &mee_fb); t->SetBranchAddress("mpipi_fb", &mpipi_fb);
        std::unordered_set<int> closed;
        int last = -2147483647;
        double grpP = 0;    // S0 rows in P of the current group (us)
        const Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            nRows[k]++;
            if (grp != last) {
                if (k == 0) { acc[0].winVarGrp += grpP * grpP; grpP = 0; }
                if (last != -2147483647) closed.insert(last);
                if (closed.count(grp)) nAdjViol[k]++;
                last = grp;
            }
            // dz_pi column = max(|pip_z - vz|, |pim_z - vz|)
            if (pip_z == pip_z && pim_z == pim_z) {
                const float d = std::max(std::fabs(pip_z - vz), std::fabs(pim_z - vz));
                nDzChk[k]++;
                if (std::fabs(d - dz_pi) > 1e-4f) nDzBad[k]++;
            }
            const bool base = isBest >= 0 && vz > -500.f;
            const bool v160 = vz > -160.f;
            const bool piv = pipv == 1 && pimv == 1;
            bool pass[5];
            pass[0] = base && cor_ok == 1 && v160 && piv;
            pass[1] = base && cor_ok == 1 && v160;
            pass[2] = base && cor_ok == 1;
            pass[3] = base && v160 && piv;                 // cor_ok fallback: every cor_ok, *_fb values for cor_ok == 0
            pass[4] = pass[0] && isBest == 1;
            for (int s = 0; s < 5; ++s) {
                if (!pass[s]) continue;
                const bool fb = (s == 3 && cor_ok != 1);
                const double M = fb ? m4_fb : m4, E = fb ? mee_fb : mee, PP = fb ? mpipi_fb : mpipi;
                if (!(M == M)) continue;
                DSelAcc& a = acc[s];
                if (M < 0.30) a.below030[k]++;
                const int fi = Fine(M);
                if (fi >= 0) a.fine[k][fi]++;
                const int sl = SliceOf(PP);
                // windows
                auto in = [&](double lo, double hi) { return M >= lo && M < hi; };
                const bool x17m = (E >= 0.0155 && E < 0.0245);
                if (in(0.4550, 0.5225)) a.win[k][wB]++;
                if (in(0.53364, 0.55416)) { a.win[k][wP]++; if (x17m) a.win[k][wX17P]++; }
                if (in(0.5312, 0.5566) && x17m) a.win[k][wX17c]++;
                if (in(0.565, 0.700)) a.win[k][wCRU]++;
                if (in(0.400, 0.5225) && sl == 4) a.win[k][wCRK]++;
                if (in(0.5225, 0.5300)) a.win[k][wCRL]++;
                if (in(0.0, 0.455)) a.win[k][wR0]++;
                if (sl == 1 && E >= 0.010 && in(0.540, 0.5425)) a.f10[k]++;
                if (in(0.5425, 0.545) && E >= 0.030 && E < 0.032) a.x2d[k]++;
                if (s == 0 && k == 0 && in(0.53364, 0.55416)) {
                    grpP += 1;
                    if (file_id == 47 || file_id == 48) a.perIIA_P++;
                    if ((fold & 1) == 0) a.fold0_P++;
                }
                if (k > 0 && fi >= 0) {
                    const int b = fi / 20;
                    a.bandAll[b]++;
                    if (sl >= 0) {
                        a.bandSlice[b]++;
                        if (sl == 4) a.bandS5[b]++;
                        if (E == E) a.bandCell[2 * sl + (E < 0.010 ? 0 : 1)][b]++;
                    }
                    if (x17m) a.bandX17[b]++;
                }
            }
        }
        if (k == 0) acc[0].winVarGrp += grpP * grpP;
        t->ResetBranchAddresses();
    }
    f->Close();
    for (int k = 0; k < 3; ++k)
        Note(Form("tree %s: rows %lld; event groups re-opened after another group (non-adjacent rows): %lld; dz_pi column vs max(|pip_z-vz|,|pim_z-vz|): %lld of %lld differ by > 1e-4 mm",
                  trees[k], nRows[k], nAdjViol[k], nDzBad[k], nDzChk[k]));
    for (int s = 0; s < 5; ++s)
        Note(Form("%-9s rows with M4 < 0.30 GeV/c2 (outside the CB grid, inside window R0 / tot): us %.0f pp %.0f mm %.0f", stag[s], acc[s].below030[0],
                  acc[s].below030[1], acc[s].below030[2]));

    // ---- comparisons
    for (int s = 0; s < 5; ++s) {
        DSelAcc& a = acc[s];
        const std::vector<double> cb = FineCB(a.fine[1], a.fine[2]);
        auto frac = [&](const std::vector<double>& num, const std::vector<double>& den) {
            std::vector<double> r(16, 0.);
            for (int b = 0; b < 16; ++b) r[b] = den[b] > 0 ? num[b] / den[b] : 0.;
            return r;
        };
        const std::vector<double> fS5 = frac(a.bandS5, a.bandSlice), fX17 = frac(a.bandX17, a.bandAll), fCell = frac(a.bandCell[3], a.bandSlice);
        TH1* hN = H(hf, Form("data/%s/all/N_win", stag[s]));
        TH1* hCB = H(hf, Form("data/%s/all/CB_win", stag[s]));
        TH1* hPP = H(hf, Form("data/%s/all/Npp_win", stag[s]));
        TH1* hMM = H(hf, Form("data/%s/all/Nmm_win", stag[s]));
        if (!hN || !hCB || !hPP || !hMM) continue;
        Note(Form("--- data %s", stag[s]));
        const int ws[8] = {wB, wP, wX17c, wCRU, wCRK, wCRL, wR0, wX17P};
        const char* wn[8] = {"B", "P", "X17c", "CRU", "CRK", "CRL", "R0", "X17P"};
        for (int q = 0; q < 8; ++q) {
            const int w = ws[q];
            Cmp(Form("%s N   %s", stag[s], wn[q]), a.win[0][w], hN->GetBinContent(w + 1), 0.5);
            Cmp(Form("%s N++ %s", stag[s], wn[q]), a.win[1][w], hPP->GetBinContent(w + 1), 0.5);
            Cmp(Form("%s N-- %s", stag[s], wn[q]), a.win[2][w], hMM->GetBinContent(w + 1), 0.5);
        }
        const double cbB = WinCB(cb, 0.4550, 0.5225), cbP = WinCB(cb, 0.53364, 0.55416), cbU = WinCB(cb, 0.565, 0.700);
        const double cbK = WinCB(cb, 0.400, 0.5225, &fS5), cbX = WinCB(cb, 0.53364, 0.55416, &fX17), cbXc = WinCB(cb, 0.5312, 0.5566, &fX17);
        const double cbR0 = WinCB(cb, 0.30, 0.455);
        Cmp(Form("%s CB  B", stag[s]), cbB, hCB->GetBinContent(wB + 1), 1e-6 * std::max(1., cbB));
        Cmp(Form("%s CB  P", stag[s]), cbP, hCB->GetBinContent(wP + 1), 1e-6 * std::max(1., cbP));
        Cmp(Form("%s CB  CRU", stag[s]), cbU, hCB->GetBinContent(wCRU + 1), 1e-6 * std::max(1., cbU));
        Cmp(Form("%s CB  CRK (x f(S5 | band))", stag[s]), cbK, hCB->GetBinContent(wCRK + 1), 1e-6 * std::max(1., cbK));
        Cmp(Form("%s CB  X17P (x f(M(ee) 15.5-24.5 | band))", stag[s]), cbX, hCB->GetBinContent(wX17P + 1), 1e-6 * std::max(1., cbX));
        Cmp(Form("%s CB  X17c", stag[s]), cbXc, hCB->GetBinContent(wX17c + 1), 1e-6 * std::max(1., cbXc));
        Cmp(Form("%s CB  R0 (grid 0.30-0.455)", stag[s]), cbR0, hCB->GetBinContent(wR0 + 1), 1e-6 * std::max(1., cbR0));
        Note(Form("%s: N - CB in B = %.1f, in P = %.1f; CB/N in P = %.5f; window geometric means 2 sqrt(N++ N--): P %.1f, X17P %.1f, CRK %.1f, B %.1f",
                  stag[s], a.win[0][wB] - cbB, a.win[0][wP] - cbP, cbP / a.win[0][wP], 2 * std::sqrt(a.win[1][wP] * a.win[2][wP]),
                  2 * std::sqrt(a.win[1][wX17P] * a.win[2][wX17P]), 2 * std::sqrt(a.win[1][wCRK] * a.win[2][wCRK]),
                  2 * std::sqrt(a.win[1][wB] * a.win[2][wB])));
        if (s == 0) {
            // F10 cell S2E2, M4 0.540-0.5425 and X2D cell (unblinded M(ee) bin 30-32 MeV)
            TH1* f10n = H(hf, "data/S0/all/N_F10_S2E2");
            TH1* f10p = H(hf, "data/S0/all/Npp_F10_S2E2");
            TH1* f10m = H(hf, "data/S0/all/Nmm_F10_S2E2");
            TH1* f10c = H(hf, "data/S0/all/CB_F10_S2E2");
            if (f10n && f10p && f10m && f10c) {
                const int b = f10n->FindBin(0.54125);
                Note(Form("F10 bin %d = [%.4f, %.4f)", b, f10n->GetXaxis()->GetBinLowEdge(b), f10n->GetXaxis()->GetBinUpEdge(b)));
                Cmp("S0 N   F10_S2E2 [0.540,0.5425)", a.f10[0], f10n->GetBinContent(b), 0.5);
                Cmp("S0 N++ F10_S2E2 [0.540,0.5425)", a.f10[1], f10p->GetBinContent(b), 0.5);
                Cmp("S0 N-- F10_S2E2 [0.540,0.5425)", a.f10[2], f10m->GetBinContent(b), 0.5);
                const double cbc = WinCB(cb, 0.540, 0.5425, &fCell);
                Cmp("S0 CB  F10_S2E2 [0.540,0.5425) (x f(S2E2 | band))", cbc, f10c->GetBinContent(b), 1e-6 * std::max(1., cbc));
            }
            TH2* x2 = (TH2*)H(hf, "data/S0/all/N_X2D");
            TH2* x2p = (TH2*)H(hf, "data/S0/all/Npp_X2D");
            if (x2 && x2p) {
                const int bx = x2->GetXaxis()->FindBin(0.54375), by = x2->GetYaxis()->FindBin(0.031);
                Cmp("S0 N   X2D [0.5425,0.545) x M(ee) [30,32) MeV", a.x2d[0], x2->GetBinContent(bx, by), 0.5);
                Cmp("S0 N++ X2D [0.5425,0.545) x M(ee) [30,32) MeV", a.x2d[1], x2p->GetBinContent(bx, by), 0.5);
                // blinding
                int nm = 0, nmBad = 0;
                for (int y = 1; y <= x2->GetNbinsY(); ++y) {
                    const double lo = x2->GetYaxis()->GetBinLowEdge(y), hi = x2->GetYaxis()->GetBinUpEdge(y);
                    const bool shouldMask = hi > 0.013 + 1e-9 && lo < 0.027 - 1e-9;
                    for (int x = 0; x <= x2->GetNbinsX() + 1; ++x) {
                        const bool masked = x2->GetBinContent(x, y) == -1.;
                        if (shouldMask) { ++nm; if (!masked) ++nmBad; }
                        else if (masked) ++nmBad;
                    }
                }
                Note(Form("blinding X2D (S0 all): %d cells must be masked (M(ee) bins overlapping 13-27 MeV, all M4 incl. flows); wrong cells: %d", nm, nmBad));
                if (nmBad) ++gNfail;
            }
            for (int r = 2; r <= 4; ++r) {
                TH1* bm = H(hf, Form("data/S0/all/N_B_mee_R%d", r));
                if (!bm) continue;
                int bad = 0;
                for (int b = 1; b <= bm->GetNbinsX(); ++b) {
                    const double lo = bm->GetXaxis()->GetBinLowEdge(b), hi = bm->GetXaxis()->GetBinUpEdge(b);
                    const bool sm = hi > 0.013 + 1e-9 && lo < 0.027 - 1e-9;
                    if (sm != (bm->GetBinContent(b) == -1.)) ++bad;
                }
                Note(Form("blinding B_mee_R%d (S0 all): wrong cells %d", r, bad));
                if (bad) ++gNfail;
            }
            // subsets and event-group error
            TH1* hIIA = H(hf, "data/S0/perIIA/N_win");
            TH1* hF0 = H(hf, "data/S0/fold0/N_win");
            if (hIIA) Cmp("S0 N P, period IIA (days 047-048)", a.perIIA_P, hIIA->GetBinContent(wP + 1), 0.5);
            if (hF0) Cmp("S0 N P, fold0 (even fold)", a.fold0_P, hF0->GetBinContent(wP + 1), 0.5);
            Cmp("S0 N P event-group error", std::sqrt(a.winVarGrp), hN->GetBinError(wP + 1), 1e-6 * std::sqrt(a.winVarGrp));
        }
    }
}

// ---------------------------------------------------------------------------- simulation
struct SQ {
    std::string name, path;
    int bin;            // bin of the file histogram (win: id + 1; tot: 1); for TH1 F10 cells: computed
    double sum = 0;
    std::unordered_map<long long, double> grp;   // event-group sums (only when errors wanted)
    bool err = false;
};
inline double MM4at(const simeng::Row& r, double T) {
    const double mp = 938.27231, me = 0.51099895, mpi = 139.56995;
    auto p4 = [](double p, double th, double ph, double m) {
        const double t = th * TMath::DegToRad(), f = ph * TMath::DegToRad();
        TLorentzVector v; v.SetXYZM(p * std::sin(t) * std::cos(f), p * std::sin(t) * std::sin(f), p * std::cos(t), m); return v;
    };
    const double E = T + mp;
    const TLorentzVector ini(0, 0, std::sqrt(E * E - mp * mp), E + mp);
    const TLorentzVector fs = p4(r.ep_p, r.ep_theta, r.ep_phi, me) + p4(r.em_p, r.em_theta, r.em_phi, me) + p4(r.pip_p, r.pip_theta, r.pip_phi, mpi) +
                              p4(r.pim_p, r.pim_theta, r.pim_phi, mpi);
    const double m2 = (ini - fs).M2();
    return (m2 >= 0 ? std::sqrt(m2) : -std::sqrt(-m2)) * 1e-3;
}

void SimPart(TFile* hf) {
    Note("=============================== SIMULATION (original FAT_simwt headers, streaming, own definitions)");
    simphys::Init();
    simphys::LoadAngles();
    const double Ls = 221.05, s = 0.99615;
    enum { qT12pre, qC23pre_s1, qC5P, qC23B, qGetaP, qT12f10, qT12modeB, qT12ghostA, qC6X17c, qT12peak, qC5P_r089, qC5P_V2, qT12wdb, qT12CRU,
           qGconvCRU, qNQ };
    std::vector<SQ> Q(qNQ);
    Q[qT12pre] = {"pre  T12 all M4 (acceptance 1: 701.1)", "sim/pre__nom/T12_tot", 1};
    Q[qC23pre_s1] = {"pre  C2+C3, sim-scale m4 0.4576-0.5255 (acc. 1: 9.689)", "", 0};
    Q[qC5P] = {"S0   C5 in P (M4')", "sim/S0__nom/C5_win", wP + 1};
    Q[qC23B] = {"S0   C2+C3 in B (M4')", "", wB + 1};
    Q[qGetaP] = {"S0   G_eta in P", "sim/S0__nom/Geta_win", wP + 1};
    Q[qT12f10] = {"S0   T12 F10_S3E1 [0.5425,0.545)", "sim/S0__nom/T12_F10_S3E1", -1};
    Q[qT12modeB] = {"S0   modeB T12 all M4", "sim/S0__modeB/T12_tot", 1};
    Q[qT12ghostA] = {"S0   ghostA T12 all M4", "sim/S0__ghostA/T12_tot", 1};
    Q[qC6X17c] = {"S0   C6 (X17, per unit BR) in X17c", "sim/S0__nom/C6_win", wX17c + 1};
    Q[qT12peak] = {"PEAK T12 all M4 (MM4 at 4500 MeV)", "sim/PEAK__nom/T12_tot", 1};
    Q[qC5P_r089] = {"S0   r = 0.89 C5 in P", "sim/S0__r089/C5_win", wP + 1};
    Q[qC5P_V2] = {"S0   V2 (107) C5 in P", "sim/S0__V2/C5_win", wP + 1};
    Q[qT12wdb] = {"S0   w_db T12 all M4", "sim/S0__wdb/T12_tot", 1};
    Q[qT12CRU] = {"S0   T12 in CR-U", "sim/S0__nom/T12_win", wCRU + 1};
    Q[qGconvCRU] = {"S0   G_conv in CR-U", "sim/S0__nom/Gconv_win", wCRU + 1};
    Q[qC5P].err = Q[qC23B].err = Q[qT12CRU].err = true;
    struct TF { const char* file; const char* tree; int tid; };
    const TF tf[3] = {{kSim, "std_us", 0}, {kSim, "lep_us", 3}, {kSimEta, "eta_us", 6}};
    long long nPeakChk = 0;
    for (const TF& x : tf) {
        TFile* f = TFile::Open(x.file, "READ");
        TTree* t = (TTree*)f->Get(x.tree);
        simeng::Row r;
        simeng::BindRow(t, r, false);
        const Long64_t n = t->GetEntries();
        long long used = 0;
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            if (!r.ownA && !r.ownB) continue;
            const bool pt3 = (r.trigbit & 8192) != 0;
            const bool pre = r.isBest >= 0 && r.vz > -500.f && r.vz > -160.f && r.cor_ok == 1 && pt3;
            if (!pre) continue;
            ++used;
            const bool piv = r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1;
            const bool sig = r.tmpl == 1 && r.cls >= 2 && r.cls <= 5;
            const bool t12 = sig || r.tmpl == 2;
            const bool gp = simeng::GhostPair(r);
            const long long gk = ((long long)x.tid << 40) | (long long)(unsigned)r.grp;
            auto add = [&](int q, double w) {
                if (w == 0) return;
                Q[q].sum += w;
                if (Q[q].err) Q[q].grp[gk] += w;
            };
            const double m4p = (r.m4_true > 0) ? s * (r.m4_true + (r.m4 - r.m4_true)) : s * r.m4;
            const double m4r = (r.m4_true > 0) ? s * (r.m4_true + 0.89 * (r.m4 - r.m4_true)) : s * r.m4;
            const double mpp = 0.9952 * r.mpipi;
            const bool inP = m4p >= 0.53364 && m4p < 0.55416, inB = m4p >= 0.4550 && m4p < 0.5225, inU = m4p >= 0.565 && m4p < 0.700;
            double wA = 0, wB_ = 0, wdb = 0, f108 = 0;
            if (r.ownA || r.ownB) f108 = simphys::RowFactor(r, 108);
            if (r.ownA) { wA = r.w * f108; wdb = r.w_db * f108; }
            if (r.ownB) wB_ = r.w * f108;
            if (r.ownA) {
                if (t12) add(qT12pre, wA);
                if ((r.tmpl == 1 && (r.cls == 2 || r.cls == 3)) && r.m4 >= 0.4576f && r.m4 < 0.5255f) add(qC23pre_s1, wA);
                // PEAK: pre + dz_pi < 40 + MM4(4500) < 2.85
                const float dz = std::max(std::fabs(r.pip_z - r.vz), std::fabs(r.pim_z - r.vz));
                if (t12 && dz < 40.f) {
                    const double mm = MM4at(r, 4500.);
                    if ((float)mm < 2.85f) add(qT12peak, wA);
                    ++nPeakChk;
                }
            }
            if (!piv) continue;
            // S0
            if (r.ownA) {
                if (r.tmpl == 1 && r.cls == 5 && inP) add(qC5P, wA);
                if (r.tmpl == 1 && (r.cls == 2 || r.cls == 3) && inB) add(qC23B, wA);
                if (r.tmpl == 2 && (r.pcls == 2 || r.lcls == 2) && inP) add(qGetaP, wA);
                if (t12 && m4p >= 0.5425 && m4p < 0.545 && mpp >= 0.37 && mpp < 0.39 && r.mee < 0.010) add(qT12f10, wA);
                if (t12 || gp) add(qT12ghostA, wA);
                if (r.tmpl == 1 && r.cls == 6 && m4p >= 0.5312 && m4p < 0.5566 && r.mee >= 0.0155 && r.mee < 0.0245) add(qC6X17c, 2.);
                if (r.tmpl == 1 && r.cls == 5) {
                    if (m4r >= 0.53364 && m4r < 0.55416) add(qC5P_r089, wA);
                    if (inP) add(qC5P_V2, r.w * simphys::RowFactor(r, 107));
                }
                if (t12) add(qT12wdb, wdb);
                if (t12 && inU) add(qT12CRU, wA);
                if (r.tmpl == 2 && !(r.pcls == 2 || r.lcls == 2) && r.lcls == 3 && inU) add(qGconvCRU, wA);
            }
            if (r.ownB && t12) add(qT12modeB, wB_);
        }
        Note(Form("sim tree %s: %lld rows, %lld owned pre rows", x.tree, (long long)n, used));
        t->ResetBranchAddresses();
        f->Close();
    }
    // comparisons
    for (int q = 0; q < qNQ; ++q) {
        SQ& x = Q[q];
        const double mine = x.sum / Ls;
        double err = 0;
        if (x.err) { for (auto& g : x.grp) err += g.second * g.second; err = std::sqrt(err) / Ls; }
        if (q == qC23pre_s1) { Cmp(x.name.c_str(), mine, 9.6887, 5e-4, "(ref: WP0 note / sim study 9.689)"); continue; }
        if (q == qC23B) {
            TH1* a = H(hf, "sim/S0__nom/C2_win");
            TH1* b = H(hf, "sim/S0__nom/C3_win");
            if (a && b) {
                Cmp(x.name.c_str(), mine, a->GetBinContent(wB + 1) + b->GetBinContent(wB + 1), 1e-9 + 1e-9 * mine);
                Cmp("S0   C2+C3 in B event-group error", err, std::sqrt(std::pow(a->GetBinError(wB + 1), 2) + std::pow(b->GetBinError(wB + 1), 2)), 0.02 * err,
                    "(C2 and C3 rows of one group combined vs separate: expect a small difference)");
            }
            continue;
        }
        TH1* h = H(hf, x.path.c_str());
        if (!h) continue;
        const int b = x.bin > 0 ? x.bin : h->FindBin(0.54375);
        Cmp(x.name.c_str(), mine, h->GetBinContent(b), 1e-9 + 1e-9 * std::fabs(mine));
        if (x.err) Cmp((x.name + " event-group error").c_str(), err, h->GetBinError(b), 1e-6 * err);
    }
    Note(Form("PEAK: MM4 at 4500 MeV computed for %lld rows (own TLorentzVector code)", nPeakChk));
}
// ---------------------------------------------------------------------------- internal consistency of the file
double Sum(TH1* h, int b1, int b2) { double s = 0; for (int b = b1; b <= b2; ++b) s += h->GetBinContent(b); return s; }
void FilePart(TFile* hf) {
    Note("=============================== FILE (internal consistency of hists_v1.root)");
    // closure F10 vs fine (data S0 N and CB, sim S0__nom T12)
    const char* blk[3] = {"data/S0/all/N_", "data/S0/all/CB_", "sim/S0__nom/T12_"};
    for (int k = 0; k < 3; ++k) {
        double s10 = 0;
        for (int sl = 1; sl <= 5; ++sl)
            for (int e = 1; e <= 2; ++e) { TH1* h = H(hf, Form("%sF10_S%dE%d", blk[k], sl, e)); if (h) s10 += Sum(h, 1, h->GetNbinsX()); }
        TH1* fi = H(hf, Form("%sfine_M4", blk[k]));
        if (fi) {
            const int b1 = fi->FindBin(0.4001), b2 = fi->FindBin(0.6999);
            Cmp(Form("closure %sF10 (10 cells) vs fine 0.40-0.70", blk[k]), s10, Sum(fi, b1, b2), 1e-9 * std::max(1., s10));
        }
    }
    // sim component sums (S0__nom, tot): T12 = C2..C5 + Geta + Gpi0D + Gconv + Gvec + Goth; GD; om/nom split
    auto tot = [&](const char* c) { TH1* h = H(hf, Form("sim/S0__nom/%s_tot", c)); return h ? h->GetBinContent(1) : 0.; };
    const double t12 = tot("T12");
    const double parts = tot("C2") + tot("C3") + tot("C4") + tot("C5") + tot("Geta") + tot("Gpi0D") + tot("Gconv") + tot("Gvec") + tot("Goth");
    Cmp("sim S0__nom T12 = C2..C5 + G_eta + G_pi0D + G_conv + G_vec + G_oth", parts, t12, 1e-9 * t12);
    Cmp("sim S0__nom GD = G_pi0D + G_vec + G_oth", tot("Gpi0D") + tot("Gvec") + tot("Goth"), tot("GD"), 1e-9 * tot("GD"));
    Cmp("sim S0__nom Gconv = _om + _nom", tot("Gconv_om") + tot("Gconv_nom"), tot("Gconv"), 1e-9 * tot("Gconv"));
    // halves
    TH1* c5 = H(hf, "sim/S0__nom/C5_win");
    TH1* c5a = H(hf, "sim/S0__nom/C5_win_hA");
    TH1* c5b = H(hf, "sim/S0__nom/C5_win_hB");
    if (c5 && c5a && c5b) {
        Cmp("sim S0__nom C5 P: half A + half B = full", c5a->GetBinContent(wP + 1) + c5b->GetBinContent(wP + 1), c5->GetBinContent(wP + 1), 1e-12);
        Note(Form("sim S0__nom C5 P: half A %.5f, half B %.5f (ratio %.4f)", c5a->GetBinContent(wP + 1), c5b->GetBinContent(wP + 1),
                  c5a->GetBinContent(wP + 1) / c5b->GetBinContent(wP + 1)));
    }
    // replicas: data CB (win P, B, CRK, X17P) and sim C5 (win P)
    TH2* rep = (TH2*)H(hf, "data/S0/all/CBrep_win");
    TH1* cbw = H(hf, "data/S0/all/CB_win");
    const int ws[4] = {wP, wB, wCRK, wX17P};
    const char* wn[4] = {"P", "B", "CRK", "X17P"};
    if (rep && cbw)
        for (int q = 0; q < 4; ++q) {
            double m = 0, m2 = 0;
            const int nr = rep->GetNbinsY();
            for (int r = 1; r <= nr; ++r) { const double x = rep->GetBinContent(ws[q] + 1, r); m += x; m2 += x * x; }
            m /= nr; m2 = std::sqrt(std::max(0., m2 / nr - m * m));
            Note(Form("data S0 CB %-4s: CB %.2f, linear error %.2f; %d replicas: mean %.2f (pull of mean %.2f), rms %.2f (rms/linear %.2f)", wn[q],
                      cbw->GetBinContent(ws[q] + 1), cbw->GetBinError(ws[q] + 1), nr, m, (m - cbw->GetBinContent(ws[q] + 1)) / (m2 / std::sqrt(nr)), m2,
                      m2 / cbw->GetBinError(ws[q] + 1)));
        }
    TH2* srep = (TH2*)H(hf, "sim/S0__nom/C5_win_rep");
    if (srep && c5) {
        double m = 0, m2 = 0;
        const int nr = srep->GetNbinsY();
        for (int r = 1; r <= nr; ++r) { const double x = srep->GetBinContent(wP + 1, r); m += x; m2 += x * x; }
        m /= nr; m2 = std::sqrt(std::max(0., m2 / nr - m * m));
        Note(Form("sim S0__nom C5 P: %.5f +- %.5f (event groups); %d replicas: mean %.5f, rms %.5f (rms/error %.3f)", c5->GetBinContent(wP + 1),
                  c5->GetBinError(wP + 1), nr, m, m2, m2 / c5->GetBinError(wP + 1)));
    }
    // unit conversion: expected raw rows of C5 in P
    if (c5) Note(Form("C5 in P, S0: %.5f nb x 6079.2 x 63/64 = %.1f expected raw rows at SMASH x V3", c5->GetBinContent(wP + 1),
                      c5->GetBinContent(wP + 1) * 6079.2 * 63. / 64.));
}
}  // namespace v0

void v0_verify(int what = 7) {
    using namespace v0;
    gOut = fopen(kOut, "w");
    fprintf(gOut, "# v0_verify.txt -- plots/all_report/v0_verify.C: independent recomputation of hists_v1.root contents (WP0 verification)\n");
    fprintf(gOut, "# data: raw TTree reading of scan_v2/scan_compact.root; simulation: original FAT_simwt headers, streaming; nb = sum w / 221.05\n");
    TFile* hf = TFile::Open(kHists, "READ");
    if (what & 1) DataPart(hf);
    if (what & 2) SimPart(hf);
    if (what & 4) FilePart(hf);
    Note(Form("TOTAL mismatches: %d", gNfail));
    fclose(gOut);
}
