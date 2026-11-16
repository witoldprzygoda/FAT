// ============================================================================
// x17ab_toys.C  --  WP9 X-b: validation of the X17 2D likelihood on the 1500 X17 toys of the simulation study plus a
//                   CB-like component (SPEC report_all)
// ============================================================================
// Inputs (read only):
//   toys     /hdd3/.../sim/pippimepem_betanew/blind/pseudo/x17_toys.root: 1500 toys at L = 6080.5 nb^-1 (the simulation
//            study's value), rows after pre + both pions in the vertex fit (= S0), simulated M4 (no mass transfer),
//            M4 0.4830-0.6080 (2.5 MeV) x M(ee) 0-100 MeV (0.5 MeV); expectation = half B (odd event groups) x 2 of
//            T1 (classes 2-5, no X17) + T2, plus X17 (class 6, m = 17 MeV/c2) at BR in {0, 1, 3, 10, 30} x 1e-6 (300 each)
//   key      .../blind/sealed/x17_key.root (injected BR per toy id; the simulation blind test was unblinded on
//            2026-11-13, tables/s5_blind_eval_x17.txt; the key is read only after all fits of a configuration are done)
//   templates .../blind/templates/templates.root (half A x 2, m4mee_X17_*): the model templates (never half B)
//   CB-like  hists_v1.root data/S0/all: CB_X2D (truth), CBrep_X2D (200 Poisson event-group replicas of the like-sign rows:
//            the "measured" CB of toy t is replica t mod 200), Npp_X2D / Nmm_X2D (geometric-mean variant)
// Toy grid: simulated M4 0.5230-0.6080 (34 bins of 2.5 MeV) x M(ee) 2-62 MeV/c2 (30 bins of 2 MeV) = the X2D layout
//   shifted to the simulation scale; the data CB of X2D bin (i, j) is used for toy bin (i, j) (offset 1.5 MeV in M4, the CB
//   is smooth in M4).
// Configurations (argument cfg):
//   0 NOMINAL: inject Poisson(CB_spread), fit with the CB replica; k_E1 1 +- 0.5, k_E2 1 +- 0.05; a_eta 1 +- 0.10; g 1 +- 0.30;
//     background B (1 + (b + b'u) t + (c + c'u) t^2) per class (E1, E2), u = (M(ee) - 30 MeV)/30 MeV, b, c, b', c' free
//     (bkgMode 4), B = exp(pol2)-smoothed half-A (G_pi0D + G_conv + G_vec + G_oth); conversion-tail
//     nuisance (M(ee) >= 12 MeV) with the half-A MC error; X = half-A C6 line
//   1 no CB at all (the toys as they are; CB term 0)
//   2 CB construction mismatch: inject Poisson(CB_geo) (per-row geometric mean 2 sqrt(N++ N--), M4 shape of CB_fine),
//     fit with the spread-rule CB replica
//   3 background "pol2": flat x (1 + b t + c t^2) per class
//   4 three background classes (E1, 10-30, >= 30 MeV)
//   5 background shape fixed (no tilt)
//   6 no conversion-tail nuisance
//   7 X line from the C5 slice at 17 MeV/c2 reweighted to X17 kinematics (raw slice of x17_lines.root), not C6
//   8 CB fixed at the nominal (no like-sign fluctuation in the fit input)
//   9 nominal, but the X17 yield nuisance sX = 0.10 (eps_X/eps_5 and N5 bias in the likelihood)
//  10 SMOOTH-PARENT toys (own, 300 per BR in {0, 1, 3, 10, 30} x 1e-6): truth = the model itself with the parameters of a
//     fit to the half-A Asimov set (sumT12 x 27.5) at r = 0, except N5, N4 = the half-A expectations and a_eta = g = 1;
//     X17 added as r N5 R_X gX X with r = BR / 2.68e-4; Poisson rows + Poisson(CB_spread); fit as cfg 0. These toys are free
//     of the heavy-row spikes of the simulated parent and test the statistics of the method (linearity, error, CLs coverage)
//  11 smooth-parent toys with the CB mismatch of cfg 2
//  12 smooth-parent toys, truth and fit with the tilt-only background B (1 + b t)
//  13 smooth-parent toys whose truth uses THREE background classes (E1, 10-30, >= 30 MeV, each its own smoothed shape),
//     fitted with the two-class model B (1 + b t + c t^2): test of a background-shape change with M(ee)
//  14 truth as 13, fit with the nominal model (two classes, bkgMode 4)
//  15 truth as 10 (two classes), fit with bkgMode 3 (no M(ee) dependence of the shape correction)
//  16 truth as 13, fit with three classes (bkgMode 3)
//  17 truth as 10, fit with three classes (bkgMode 3)
//  18 1500 simulation-study toys + CB-like, fit with bkgMode 3
//  3 = pol2 (bkgMode 2), 4 = three classes with bkgMode 3, 5 = fixed shape (bkgMode 1)
// Nominal k_E1: 1 +- 0.5 (see the results note: free k_E1 is degenerate with the free n_j of the E1 rows and ran away to -19)
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/x17ab/toys_cfg<cfg>.root (tree "fits", one entry per toy)
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17ab", kTRUE)' \
//        'plots/all_report/x17ab_toys.C+O(0, 0, 1500)'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"

#include "defs.h"
#include "hists_io.h"
#include "x17_model.h"

namespace x17toys {
const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/x17ab";

// data CB (and like-sign) of X2D mapped onto the toy grid
struct CBset {
    std::vector<double> nom, geo;
    std::vector<std::vector<double>> rep;
};
bool LoadCB(const x17m::Inputs& I, CBset& S) {
    TFile* f = TFile::Open(alldefs::kHistFile, "READ");
    if (!f || f->IsZombie()) return false;
    TH2* cb = (TH2*)f->Get("data/S0/all/CB_X2D");
    TH3* rp = (TH3*)f->Get("data/S0/all/CBrep_X2D");
    TH2* pp = (TH2*)f->Get("data/S0/all/Npp_X2D");
    TH2* mm = (TH2*)f->Get("data/S0/all/Nmm_X2D");
    if (!cb || !rp || !pp || !mm) { printf("x17toys: CB histograms missing\n"); return false; }
    const int nx = I.G.nx, ny = I.G.ny;
    S.nom.assign(I.G.n(), 0.); S.geo.assign(I.G.n(), 0.);
    // geometric-mean variant: row total 2 sqrt(N++_j N--_j) over the X2D M4 range, M4 shape of CB_fine (sum over M(ee) incl. flows)
    std::vector<double> fine(cb->GetNbinsX(), 0.);
    double ftot = 0;
    for (int i = 0; i < cb->GetNbinsX(); ++i) { for (int j = 0; j <= cb->GetNbinsY() + 1; ++j) fine[i] += cb->GetBinContent(i + 1, j); ftot += fine[i]; }
    for (int j = 0; j < ny; ++j) {
        double sp = 0, sm = 0;
        for (int i = 1; i <= cb->GetNbinsX(); ++i) { sp += pp->GetBinContent(i, j + 1); sm += mm->GetBinContent(i, j + 1); }
        const double row = 2 * std::sqrt(sp * sm);
        for (int i = 0; i < nx; ++i) {
            S.nom[I.G.idx(i, j)] = cb->GetBinContent(i + 1, j + 1);
            S.geo[I.G.idx(i, j)] = row * fine[i] / ftot;
        }
    }
    const int nr = rp->GetNbinsZ();
    S.rep.assign(nr, std::vector<double>(I.G.n(), 0.));
    for (int r = 0; r < nr; ++r) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) S.rep[r][I.G.idx(i, j)] = rp->GetBinContent(i + 1, j + 1, r + 1);
    f->Close(); delete f;
    return true;
}
}  // namespace x17toys

void x17ab_toys(int cfg = 0, int first = 0, int ntoy = 1500, unsigned long long seed = 0x17ab0000ULL) {
    using namespace x17toys;
    gSystem->mkdir(kOut, kTRUE);
    TStopwatch sw;
    // ---------------------------------------------------------------- model inputs on the toy grid (half A)
    x17m::ToyGridOpt to;
    const bool truth3 = (cfg == 13 || cfg == 14 || cfg == 16);
    to.nbcl = (cfg == 4 || cfg == 16 || cfg == 17) ? 3 : 2;
    x17m::Inputs I;
    if (!x17m::BuildToyGrid(I, to)) { printf("x17toys: cannot build the toy grid\n"); return; }
    CBset CB;
    if (!LoadCB(I, CB)) return;
    // conversion-tail uncertainty from the half-A template errors (M(ee) >= 12 MeV, M4 range of the grid)
    double sC = 0;
    {
        TFile* f = TFile::Open(x17m::kToyTemplates, "READ");
        TH2* c4 = (TH2*)f->Get("m4mee_X17_conv4");
        std::vector<double> e2;
        std::vector<double> v = x17m::ToyRebin(c4, to, &e2);
        double s = 0, ee = 0;
        for (int j = 0; j < I.G.ny; ++j) if (I.Ctail[j] > 0) for (int i = 0; i < I.G.nx; ++i) { s += v[I.G.idx(i, j)]; ee += e2[I.G.idx(i, j)]; }
        // the M(ee) flows above 62 MeV are part of the tail as well, but carry no conversions in practice
        sC = s > 0 ? std::sqrt(ee) / s : 0.;
        f->Close(); delete f;
    }
    if (cfg == 7) {
        std::vector<double> X; double RX;
        if (!x17m::LineFor(17.0, I.G, X, RX, 0, 1., 0., true)) { printf("x17toys: no line file\n"); return; }
        double s0 = 0, s1 = 0; for (double x : I.X) s0 += x; for (double x : X) s1 += x;
        printf("x17toys: cfg 7: X from the C5 slice at 17 (sum in 2-62 MeV %.4f; half-A C6 %.4f); R_X kept from half A\n", s1, s0);
        I.X = X;
    }
    x17m::Config C;
    C.kSig[0] = 0.5; C.kSig[1] = 0.05;
    C.a0 = 1.; C.aSig = 0.10; C.gSig = 0.30;
    // nominal background model: bkgMode 4 = B (1 + (b + b' u) t + (c + c' u) t^2), two classes
    C.bkgMode = (cfg == 3) ? 2 : (cfg == 5 ? 1 : (cfg == 12 ? 0 : ((cfg == 4 || cfg == 13 || cfg == 15 || cfg == 16 || cfg == 17 || cfg == 18) ? 3 : 4)));
    C.cTailSig = (cfg == 6) ? 0. : sC;
    C.xSig = (cfg == 9) ? 0.10 : 0.;
    if (cfg == 1) { C.kSig[0] = 0.; C.kSig[1] = 0.; }
    const bool smoothToys = cfg >= 10 && cfg != 18;   // no CB in the toys: the CB term is switched off (CB input zero, k fixed)
    printf("x17toys: cfg %d; grid %d x %d; R_X(half A) = %.4f; N5 start %.1f, N4 start %.1f; conversion-tail sigma %.3f; bkgMode %d, classes %d\n",
           cfg, I.G.nx, I.G.ny, I.RX, I.N5start, I.N4start, sC, C.bkgMode, I.nbcl);
    double sumD = 0, sumC = 0, sumX = 0;
    for (int j = 0; j < I.G.ny; ++j) { sumD += I.D[j]; sumC += I.C[j]; sumX += I.X[j]; }
    printf("x17toys: fractions in 2-62 MeV: D %.4f, C %.4f, X %.4f\n", sumD, sumC, sumX);

    // ---------------------------------------------------------------- truth of the smooth-parent toys
    std::vector<double> truthMu, truthX;
    if (smoothToys) {
        TFile* f = TFile::Open(x17m::kToyTemplates, "READ");
        std::vector<double> nA = x17m::ToyRebin((TH2*)f->Get("m4mee_X17_sumT12"), to);
        for (auto& v : nA) v *= to.scale;
        f->Close(); delete f;
        x17m::Inputs It;
        {   // truth: two or three background classes, bkgMode 3 (cfg 12: 0)
            x17m::ToyGridOpt tt = to; tt.nbcl = truth3 ? 3 : 2;
            x17m::BuildToyGrid(It, tt);
        }
        It.CB = CB.nom;
        for (size_t b = 0; b < nA.size(); ++b) nA[b] += CB.nom[b];
        x17m::Config Ct = C; Ct.kSig[0] = 0.; Ct.kSig[1] = 0.; Ct.bkgMode = (cfg == 12) ? 0 : 3;
        x17m::Model Mt(It, Ct);
        x17m::FitResult ft0 = Mt.Fit(nA, 0., false, x17m::GlobOf(Ct));
        std::vector<double> p = ft0.p;
        p[Mt.L.iN5] = I.N5start; p[Mt.L.iN4] = I.N4start; p[Mt.L.iA] = 1.; p[Mt.L.iG] = 1.; p[Mt.L.iTC] = 0.; p[Mt.L.iTX] = 0.;
        p[Mt.L.iR] = 0.;
        truthMu = Mt.Expected(p);
        for (size_t b = 0; b < truthMu.size(); ++b) truthMu[b] = std::max(0., truthMu[b] - CB.nom[b]);   // CB injected separately
        std::vector<double> p1 = p; p1[Mt.L.iR] = 1.;
        std::vector<double> m1 = Mt.Expected(p1);
        truthX.resize(m1.size());
        for (size_t b = 0; b < m1.size(); ++b) truthX[b] = m1[b] - CB.nom[b] - truthMu[b];   // X17 rows at r = 1
        double sx = 0, sm = 0; for (double v : truthX) sx += v; for (double v : truthMu) sm += v;
        printf("x17toys: smooth parent from the half-A Asimov fit (ok %d, F %.1f): %.0f rows (no CB); X17 rows per unit r %.1f (BR 1e-5: %.1f)\n",
               ft0.ok, ft0.F, sm, sx, sx * 1e-5 / x17m::kBRdir);
    }
    const double kBRset[5] = {0., 1e-6, 3e-6, 1e-5, 3e-5};

    // ---------------------------------------------------------------- toys
    TFile* ft = TFile::Open(x17m::kToyFile, "READ");
    TTree* tl = (TTree*)ft->Get("toys");
    int tid = 0, tidx = 0;
    tl->SetBranchAddress("id", &tid); tl->SetBranchAddress("index", &tidx);
    const int nAll = smoothToys ? 1500 : (int)tl->GetEntries();
    const int last = std::min(nAll, first + ntoy);
    TFile* fo = TFile::Open(Form("%s/toys_cfg%d_%d.root", kOut, cfg, first), "RECREATE");
    TTree* out = new TTree("fits", Form("X17 toy fits, configuration %d", cfg));
    int oid, ook, oit;
    double orh, orhe, oul, oulPL, osA, oN5, oN4, ok1, ok2, oa, og, otc, otx, oF0, oFh, ochi2, oNcb, oNtoy;
    int ondf;
    out->Branch("id", &oid, "id/I"); out->Branch("ok", &ook, "ok/I"); out->Branch("iter", &oit, "iter/I");
    out->Branch("rhat", &orh, "rhat/D"); out->Branch("rhatErr", &orhe, "rhatErr/D"); out->Branch("ul", &oul, "ul/D");
    out->Branch("ulPL", &oulPL, "ulPL/D"); out->Branch("sigmaA", &osA, "sigmaA/D");
    out->Branch("N5", &oN5, "N5/D"); out->Branch("N4", &oN4, "N4/D"); out->Branch("kE1", &ok1, "kE1/D"); out->Branch("kE2", &ok2, "kE2/D");
    out->Branch("aEta", &oa, "aEta/D"); out->Branch("g", &og, "g/D"); out->Branch("thC", &otc, "thC/D"); out->Branch("thX", &otx, "thX/D");
    out->Branch("F0", &oF0, "F0/D"); out->Branch("Fhat", &oFh, "Fhat/D"); out->Branch("chi2", &ochi2, "chi2/D"); out->Branch("ndf", &ondf, "ndf/I");
    out->Branch("Ncb", &oNcb, "Ncb/D"); out->Branch("Ntoy", &oNtoy, "Ntoy/D");
    double obr = -1;
    out->Branch("brTrue", &obr, "brTrue/D");   // smooth-parent toys only (-1: join the sealed key)
    std::mt19937_64 rng(seed + 1000003ULL * cfg);
    const double scaleR = 1e-5 / x17m::kBRdir;   // typical r (BR 1e-5)
    for (int t = first; t < last; ++t) {
        TH2* h = nullptr;
        std::vector<double> n;
        if (!smoothToys) {
            tl->GetEntry(t);
            h = (TH2*)ft->Get(Form("toy_%d_m4mee", tid));
            if (!h) { printf("x17toys: toy %d missing\n", tid); continue; }
            n = x17m::ToyRebin(h, to);
        } else {
            tid = 100000 + t;
            obr = kBRset[t / 300];
            n.resize(truthMu.size());
            std::mt19937_64 r1(seed * 31 + 7919ULL * (unsigned long long)t + 17ULL * cfg);
            for (size_t b = 0; b < n.size(); ++b) {
                const double mu = truthMu[b] + obr / x17m::kBRdir * truthX[b];
                if (mu <= 0) { n[b] = 0; continue; }
                std::poisson_distribution<long long> P(mu);
                n[b] = (double)P(r1);
            }
        }
        double ntoyRows = 0; for (double v : n) ntoyRows += v;
        // CB-like component
        double ncb = 0;
        x17m::Inputs It = I;
        if (cfg != 1) {
            const std::vector<double>& truth = (cfg == 2 || cfg == 11) ? CB.geo : CB.nom;
            std::mt19937_64 r2(seed ^ (0x9e3779b97f4a7c15ULL * (unsigned long long)(tid + 1)));   // same CB rows for every cfg using the same truth
            for (size_t b = 0; b < n.size(); ++b) {
                if (truth[b] <= 0) continue;
                std::poisson_distribution<long long> P(truth[b]);
                const double k = (double)P(r2);
                n[b] += k; ncb += k;
            }
            It.CB = (cfg == 8) ? CB.nom : CB.rep[t % CB.rep.size()];
        }
        x17m::Model M(It, C);
        const x17m::Glob g = x17m::GlobOf(C);
        x17m::Limit Lm = x17m::LimitAsym(M, n, g, scaleR, false);
        const x17m::FitResult& fh = Lm.prof.hat;
        oid = tid; ook = (Lm.prof.hat.ok && Lm.prof.zero.ok) ? 1 : 0; oit = fh.iter;
        orh = Lm.rhat; orhe = Lm.rhatErr; oul = Lm.ul; oulPL = Lm.ulPL; osA = Lm.sigmaA;
        const x17m::Layout& L = M.L;
        oN5 = fh.p[L.iN5]; oN4 = fh.p[L.iN4]; ok1 = fh.p[L.iK]; ok2 = fh.p[L.iK + 1]; oa = fh.p[L.iA]; og = fh.p[L.iG];
        otc = fh.p[L.iTC]; otx = fh.p[L.iTX];
        oF0 = Lm.prof.zero.F; oFh = fh.F;
        // Pearson chi2 of the free fit
        ochi2 = 0; ondf = 0;
        for (int j = 0; j < It.G.ny; ++j) for (int i = 0; i < It.G.nx; ++i) {
            const double mu = M.Mu(fh.p, i, j);
            if (mu > 0) { ochi2 += (n[It.G.idx(i, j)] - mu) * (n[It.G.idx(i, j)] - mu) / mu; ++ondf; }
        }
        int nfree = 0; { std::vector<char> fr = M.FreeMask(true); for (char c : fr) nfree += c; }
        ondf -= nfree;
        oNcb = ncb; oNtoy = ntoyRows;
        out->Fill();
        if (h) delete h;
        if ((t - first) % 50 == 0)
            printf("x17toys: cfg %d toy %d/%d id %d: rhat %.4e +- %.4e (BR %.2e), UL %.4e (BR %.2e), N5 %.0f, kE1 %.3f kE2 %.3f, chi2/ndf %.1f/%d, %.1f s\n", cfg,
                   t - first, last - first, tid, orh, orhe, orh * x17m::kBRdir, oul, oul * x17m::kBRdir, oN5, ok1, ok2, ochi2, ondf, sw.RealTime());
        sw.Continue();
    }
    out->Write();
    fo->Close();
    ft->Close();
    printf("x17toys: cfg %d done (%d toys) in %.0f s\n", cfg, last - first, sw.RealTime());
}
