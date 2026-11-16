// ============================================================================
// wp4_fill.C  --  WP4 stress closure: simulation templates on the refined toy cells
// ============================================================================
// SPEC report_all, WP4 (prefix wp4). One pass over the simulated e+e- rows
// (std_us + lep_us + eta_us, owned rows, mode A, w x V3 = physics code 108,
// selection S0), filling for each half (A = even event group, B = odd) and for
// four mass-transfer variants
//     nom  : s = 0.99615, r = 1       (nominal, = hists_v1 S0__nom)
//     f    : s = 0.99705, r = 0.89    (WP4 scenario (f): s + 0.0009 and r = 0.89 together)
//     sUp  : s = 0.99705, r = 1       (check against hists_v1 S0__sUp)
//     r089 : s = 0.99615, r = 0.89    (check against hists_v1 S0__r089)
// the components C2, C3, C4, C5, G_eta, G_D (omega/rho0 pions, others), G_C = G_conv
// (omega/rho0 pions, others) in
//     R   : refined M4 cells (wp4_grid.h) x 15 (slice x M(ee) class e3)
//     F15 : F10 M4 bins x 15                 (event-group variance per class e3)
//     F10 : F10 M4 bins x 10 (slice x E1/E2) (event-group variance as in hists_v1)
//     tot : all M4, all slices
// Units: nb per bin (sum w / L_sim), errors = event-group errors.
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/wp4/wp4_templates.root
// Checks against hists_v1.root (S0__nom halves, S0__sUp, S0__r089): F10, F5g, win, tot.
//
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp4", kTRUE)' \
//        'plots/all_report/wp4_fill.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "wp4_grid.h"
#include "sim_eta_physics.h"
#include "hists_io.h"

namespace wp4fill {
using namespace wp4;

Grid gG;
std::vector<int> gF10x;   // fine bin -> F10 bin
int nRc, nFc;             // cells per block
// block layout: R (nR+2) x 15, F15 (nF+2) x 15, F10 (nF+2) x 10, tot 1
int oR, oF15, oF10, oTot, stride;
inline size_t Base(int v, int h, int c) { return ((size_t)(v * 2 + h) * kNTComp + c) * stride; }

struct Acc {
    std::vector<double> s, v, acc;
    std::vector<int> stamp, touched;
    int gid = 0;
    void Init(size_t n) { s.assign(n, 0.); v.assign(n, 0.); acc.assign(n, 0.); stamp.assign(n, 0); touched.reserve(1 << 12); }
    inline void Begin() { Flush(); ++gid; }
    inline void Add(size_t i, double w) {
        s[i] += w;
        if (stamp[i] != gid) { stamp[i] = gid; acc[i] = 0.; touched.push_back((int)i); }
        acc[i] += w;
    }
    void Flush() { for (int i : touched) v[i] += acc[i] * acc[i]; touched.clear(); }
};

inline int MapComp(int c, int* out) {
    int n = 0;
    switch (c) {
        case kC2: out[n++] = kTC2; break;
        case kC3: out[n++] = kTC3; break;
        case kC4: out[n++] = kTC4; break;
        case kC5: out[n++] = kTC5; break;
        case kGeta: out[n++] = kTGeta; break;
        case kGpi0Dom: case kGvecom: case kGothom: out[n++] = kTGDom; break;
        case kGpi0Dnom: case kGvecnom: case kGothnom: out[n++] = kTGDnom; break;
        case kGconvom: out[n++] = kTGCom; break;
        case kGconvnom: out[n++] = kTGCnom; break;
        default: break;
    }
    return n;
}

struct Chk { std::string what; double a, b; double tol; bool ok; };
std::vector<Chk> gChk;
void AddChk(const std::string& w, double a, double b, double tol = 1e-9) {
    const double d = std::fabs(a - b) / std::max(1e-300, std::max(std::fabs(a), std::fabs(b)));
    const bool ok = (std::fabs(a) < 1e-15 && std::fabs(b) < 1e-15) || d <= tol;
    gChk.push_back({w, a, b, tol, ok});
}

}  // namespace wp4fill

void wp4_fill() {
    using namespace wp4fill;
    TStopwatch sw;
    gG.Init();
    gF10x.assign(kNfine, 0);
    for (int i = 0; i < kNfine; ++i) gF10x[i] = FindBinE(gG.f10, kFineLo + (i + 0.5) * kFineW);
    oR = 0; oF15 = oR + (gG.nR + 2) * kNh15; oF10 = oF15 + (gG.nF + 2) * kNh15; oTot = oF10 + (gG.nF + 2) * 10; stride = oTot + 1;
    printf("wp4_fill: refined bins %d, F10 bins %d, stride %d\n", gG.nR, gG.nF, stride);
    Acc A;
    A.Init((size_t)kNMVar * 2 * kNTComp * stride);

    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    simeng::SimData* D = new simeng::SimData();
    if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
    printf("wp4_fill: loaded in %.1f s\n", sw.RealTime()); sw.Continue();

    const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};
    long long nPass = 0;
    for (int t = 0; t < 3; ++t) {
        const std::vector<simeng::Row>& vr = D->t[trees[t]];
        int cur = -1, h = 0;
        for (size_t i = 0; i < vr.size(); ++i) {
            const simeng::Row& r = vr[i];
            if (r.grp != cur) { A.Begin(); cur = r.grp; h = (r.grp % 2 == 0) ? 0 : 1; }
            if (!(r.ownA == 1)) continue;
            const float dzpi = DzPi(r.pip_z, r.pim_z, r.vz);
            if (!PassSel(kSelS0, r, r.cor_ok == 1, dzpi, (float)r.mm4)) continue;
            int cl[6];
            const int nc = CompsOf(r, false, cl);
            int tc[8], ntc = 0;
            for (int q = 0; q < nc; ++q) ntc += MapComp(cl[q], tc + ntc);
            if (ntc == 0) continue;
            const double w = r.w * simphys::RowFactor(r, kPhysNom);
            if (w == 0.) continue;
            ++nPass;
            const int j = Slice(kSPiPi * r.mpipi), e3 = E3class(r.mee), E = Eclass(r.mee);
            for (int v = 0; v < kNMVar; ++v) {
                const double m4 = M4sim(r.m4, r.m4_true, MVarS(v), MVarR(v));
                const int rb = FindBinE(gG.ref, m4);
                const int fi = FineBin(m4);
                const int fb = fi >= 0 ? gF10x[fi] : (m4 < kFineLo ? 0 : gG.nF + 1);
                for (int q = 0; q < ntc; ++q) {
                    const size_t b = Base(v, h, tc[q]);
                    if (j >= 0 && e3 >= 0) {
                        A.Add(b + oR + rb + (gG.nR + 2) * H15(j, e3), w);
                        A.Add(b + oF15 + fb + (gG.nF + 2) * H15(j, e3), w);
                        A.Add(b + oF10 + fb + (gG.nF + 2) * (2 * j + E), w);
                    }
                    A.Add(b + oTot, w);
                }
            }
        }
        A.Flush();
        printf("wp4_fill: tree %s done, %.1f s\n", simeng::TreeName(trees[t]), sw.RealTime()); sw.Continue();
    }
    delete D;
    printf("wp4_fill: %lld S0 rows with a component\n", nPass);

    // ---------------------------------------------------------------- write
    gSystem->mkdir(kOutDir, kTRUE);
    TFile* fo = TFile::Open(kTmplFile, "RECREATE");
    const double sc = 1.0 / kLsim;
    auto mk = [&](const char* name, const std::vector<double>& e, size_t off, int n) {
        TH1D* hh = new TH1D(name, "", (int)e.size() - 1, e.data());
        hh->SetDirectory(nullptr); hh->Sumw2();
        for (int k = 0; k < n; ++k) { hh->SetBinContent(k, A.s[off + k] * sc); hh->SetBinError(k, std::sqrt(A.v[off + k]) * sc); }
        return hh;
    };
    const std::vector<double> e1{0., 1.};
    for (int v = 0; v < kNMVar; ++v)
        for (int h = 0; h < 2; ++h) {
            TDirectory* d = fo->mkdir(Form("%s_h%s", MVarTag(v), h == 0 ? "A" : "B"));
            for (int c = 0; c < kNTComp; ++c) {
                const size_t b = Base(v, h, c);
                for (int j = 0; j < kNj; ++j)
                    for (int e3 = 0; e3 < kNe3; ++e3) {
                        TH1D* x = mk(Form("%s_R_S%de%d", TCompTag(c), j + 1, e3), gG.ref, b + oR + (gG.nR + 2) * H15(j, e3), gG.nR + 2);
                        d->WriteTObject(x); delete x;
                        x = mk(Form("%s_F15_S%de%d", TCompTag(c), j + 1, e3), gG.f10, b + oF15 + (gG.nF + 2) * H15(j, e3), gG.nF + 2);
                        d->WriteTObject(x); delete x;
                    }
                for (int k = 0; k < 10; ++k) {
                    TH1D* x = mk(Form("%s_F10_S%dE%d", TCompTag(c), k / 2 + 1, k % 2 + 1), gG.f10, b + oF10 + (gG.nF + 2) * k, gG.nF + 2);
                    d->WriteTObject(x); delete x;
                }
                TH1D* x = new TH1D(Form("%s_tot", TCompTag(c)), "", 1, 0., 1.);
                x->SetDirectory(nullptr);
                x->SetBinContent(1, A.s[b + oTot] * sc); x->SetBinError(1, std::sqrt(A.v[b + oTot]) * sc);
                d->WriteTObject(x); delete x;
            }
        }

    // ---------------------------------------------------------------- checks against hists_v1
    TFile* fh = TFile::Open(kHistFile);
    struct Map { int comp; std::vector<int> mine; bool err; };
    const Map maps[7] = {{kC2, {kTC2}, true}, {kC3, {kTC3}, true}, {kC4, {kTC4}, true}, {kC5, {kTC5}, true},
                         {kGeta, {kTGeta}, true}, {kGD, {kTGDom, kTGDnom}, false}, {kGconv, {kTGCom, kTGCnom}, false}};
    struct VSet { int v; const char* var; const char* suf; int hmask; };   // hmask 1 = A, 2 = B, 3 = A + B
    const VSet vs[4] = {{kMnom, "nom", "_hA", 1}, {kMnom, "nom", "_hB", 2}, {kMsUp, "sUp", "", 3}, {kMr089, "r089", "", 3}};
    double worst = 0;
    for (const auto& V : vs)
        for (const auto& M : maps) {
            auto sumMine = [&](size_t off, double* var) {
                double s = 0, q = 0;
                for (int h = 0; h < 2; ++h) {
                    if (!(V.hmask & (1 << h))) continue;
                    for (int c : M.mine) { s += A.s[Base(V.v, h, c) + off]; q += A.v[Base(V.v, h, c) + off]; }
                }
                if (var) *var = q;
                return s * sc;
            };
            // F10
            for (int k = 0; k < 10; ++k) {
                TH1* hf = allio::Sim(fh, "S0", V.var, M.comp, allio::F10(k / 2, k % 2), V.suf);
                if (!hf) continue;
                for (int x = 1; x <= gG.nF; ++x) {
                    double var;
                    const double mine = sumMine(oF10 + (gG.nF + 2) * k + x, &var);
                    AddChk(Form("%s%s %s F10 k%d x%d", V.var, V.suf, CompTag(M.comp), k, x), mine, hf->GetBinContent(x));
                    worst = std::max(worst, std::fabs(mine - hf->GetBinContent(x)) / std::max(1e-12, std::fabs(hf->GetBinContent(x))));
                    if (M.err && (V.hmask != 3))
                        AddChk(Form("%s%s %s F10 err k%d x%d", V.var, V.suf, CompTag(M.comp), k, x), std::sqrt(var) * sc, hf->GetBinError(x));
                }
                // refined sum per F10 bin = F10 (content)
                delete hf;
            }
            // F5g = F15 class 2
            for (int j = 0; j < kNj; ++j) {
                TH1* hf = allio::Sim(fh, "S0", V.var, M.comp, allio::F5g(j), V.suf);
                if (!hf) continue;
                for (int x = 1; x <= gG.nF; ++x)
                    AddChk(Form("%s%s %s F5g S%d x%d", V.var, V.suf, CompTag(M.comp), j, x),
                           sumMine(oF15 + (gG.nF + 2) * H15(j, 2) + x, nullptr), hf->GetBinContent(x));
                delete hf;
            }
            // windows from refined cells (all slices and classes)
            TH1* hw = allio::Sim(fh, "S0", V.var, M.comp, "win", V.suf);
            const int wins[6] = {kWB, kWP, kWSBlo, kWSBhi, kWCRU, kWCRK};
            for (int w : wins) {
                double s = 0;
                for (int r : gG.RefIn(kWin[w].lo, kWin[w].hi))
                    for (int h15 = 0; h15 < kNh15; ++h15) {
                        if (w == kWCRK && h15 / kNe3 != 4) continue;
                        s += sumMine(oR + (gG.nR + 2) * h15 + r, nullptr);
                    }
                AddChk(Form("%s%s %s win %s", V.var, V.suf, CompTag(M.comp), kWin[w].tag), s, allio::Win(hw, w));
            }
            delete hw;
            TH1* ht = allio::Sim(fh, "S0", V.var, M.comp, "tot", V.suf);
            AddChk(Form("%s%s %s tot", V.var, V.suf, CompTag(M.comp)), sumMine(oTot, nullptr), ht->GetBinContent(1));
            delete ht;
        }
    // refined -> F10 closure (all variants, both halves, all components)
    for (int v = 0; v < kNMVar; ++v)
        for (int h = 0; h < 2; ++h)
            for (int c = 0; c < kNTComp; ++c) {
                const size_t b = Base(v, h, c);
                double maxd = 0;
                for (int h15 = 0; h15 < kNh15; ++h15) {
                    std::vector<double> s(gG.nF + 2, 0.);
                    for (int r = 1; r <= gG.nR; ++r) s[gG.refToF10[r]] += A.s[b + oR + (gG.nR + 2) * h15 + r];
                    for (int x = 1; x <= gG.nF; ++x) {
                        const double f = A.s[b + oF15 + (gG.nF + 2) * h15 + x];
                        maxd = std::max(maxd, std::fabs(s[x] - f) / std::max(1e-12, std::fabs(f)));
                    }
                }
                AddChk(Form("closure R->F15 %s h%d %s", MVarTag(v), h, TCompTag(c)), 1.0 + maxd, 1.0, 1e-9);
            }
    int nBad = 0;
    FILE* tab = fopen(Form("%s/wp4_fill_checks.txt", kOutDir), "w");
    fprintf(tab, "# wp4_fill.C checks against hists_v1.root (tolerance 1e-9 relative)\n");
    for (const auto& c : gChk) {
        if (!c.ok) ++nBad;
        if (!c.ok) fprintf(tab, "BAD  %-50s mine %.10g file %.10g\n", c.what.c_str(), c.a, c.b);
    }
    fprintf(tab, "checks %zu, bad %d, worst F10 content relative difference %.3g\n", gChk.size(), nBad, worst);
    fclose(tab);
    printf("wp4_fill: checks %zu, bad %d, worst F10 relative difference %.3g\n", gChk.size(), nBad, worst);
    int nPr = 0;
    for (const auto& c : gChk) if (!c.ok && nPr++ < 20) printf("  BAD %s mine %.10g file %.10g\n", c.what.c_str(), c.a, c.b);
    fo->Close();
    printf("wp4_fill: wrote %s in %.1f s\n", kTmplFile, sw.RealTime());
}
