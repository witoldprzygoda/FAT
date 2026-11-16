// ============================================================================
// x17ab_toycls.C  --  WP9 X-b: asymptotic against toy-based CLs for the X17 2D likelihood (SPEC report_all)
// ============================================================================
// For a few data sets the 90 % CLs limit on r is computed twice:
//   asymptotic (x17m::LimitAsym: q~_r, sigma from the background-only Asimov set), and
//   with toys (x17m::ToyCLs: for every r on a grid, Poisson pseudo-data from the conditional fits at r and at 0 with the
//   centres of the Gaussian constraints drawn around the conditional estimates; CLs = P(q~ >= q_obs | r) / P(q~ >= q_obs | 0);
//   the limit is the crossing of 0.10, log-linear interpolation).
// Data sets (argument which):
//   0 toy grid, smooth parent (as cfg 10 of x17ab_toys.C) + CB, background-only Asimov set (the median expectation)
//   1, 2 smooth-parent toys at BR = 0;  3, 4 smooth-parent toys at BR = 1e-5
//   5, 6 simulation-study toys at BR = 0 (+ CB-like), fitted with the nominal model
//   7 data-scale X2D model: the pre-data projection Asimov set of x17ab_model.C (a_eta 0.40, background 1.35 x SMASH, data CB)
// The rule of the spec: if the two limits differ by more than 10 %, the toys are primary.
// Output: /hdd3/.../x17ab/toycls_<which>.txt (one line per mu and a summary line)
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17ab", kTRUE)' 'plots/all_report/x17ab_toycls.C+O(0, 1000)'
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

namespace x17tc {
const char* const kOut = "/hdd3/przygoda/HADES/FAT_out/all_report/x17ab";
}

void x17ab_toycls(int which = 0, int ntoys = 1000) {
    using namespace x17tc;
    TStopwatch sw;
    x17m::Inputs I;
    x17m::Config C;
    C.kSig[0] = 0.5; C.kSig[1] = 0.05; C.aSig = 0.10; C.gSig = 0.30; C.bkgMode = 4;
    std::vector<double> nobs;
    std::string label;
    if (which <= 6) {
        // ---------------------------------------------------------------- toy grid
        x17m::ToyGridOpt to;
        if (!x17m::BuildToyGrid(I, to)) return;
        TFile* f = TFile::Open(alldefs::kHistFile, "READ");
        TH2* cb = (TH2*)f->Get("data/S0/all/CB_X2D");
        TH3* rp = (TH3*)f->Get("data/S0/all/CBrep_X2D");
        std::vector<double> cbn(I.G.n()), cbr(I.G.n());
        for (int j = 0; j < I.G.ny; ++j) for (int i = 0; i < I.G.nx; ++i) {
            cbn[I.G.idx(i, j)] = cb->GetBinContent(i + 1, j + 1);
            cbr[I.G.idx(i, j)] = rp->GetBinContent(i + 1, j + 1, which + 1);
        }
        f->Close(); delete f;
        TFile* ft = TFile::Open(x17m::kToyTemplates, "READ");
        TH2* c4 = (TH2*)ft->Get("m4mee_X17_conv4");
        std::vector<double> e2;
        std::vector<double> v = x17m::ToyRebin(c4, to, &e2);
        double s = 0, ee = 0;
        for (int j = 0; j < I.G.ny; ++j) if (I.Ctail[j] > 0) for (int i = 0; i < I.G.nx; ++i) { s += v[I.G.idx(i, j)]; ee += e2[I.G.idx(i, j)]; }
        C.cTailSig = s > 0 ? std::sqrt(ee) / s : 0.;
        std::vector<double> nA = x17m::ToyRebin((TH2*)ft->Get("m4mee_X17_sumT12"), to);
        for (auto& q : nA) q *= to.scale;
        ft->Close(); delete ft;
        if (which <= 4) {
            // smooth parent (as x17ab_toys.C cfg 10)
            x17m::Inputs It = I; It.CB = cbn;
            for (size_t b = 0; b < nA.size(); ++b) nA[b] += cbn[b];
            x17m::Config Ct = C; Ct.kSig[0] = 0.; Ct.kSig[1] = 0.; Ct.bkgMode = 3;
            x17m::Model Mt(It, Ct);
            x17m::FitResult f0 = Mt.Fit(nA, 0., false, x17m::GlobOf(Ct));
            std::vector<double> p = f0.p;
            p[Mt.L.iN5] = I.N5start; p[Mt.L.iN4] = I.N4start; p[Mt.L.iA] = 1.; p[Mt.L.iG] = 1.; p[Mt.L.iTC] = 0.; p[Mt.L.iR] = 0.;
            std::vector<double> mu0 = Mt.Expected(p);
            std::vector<double> p1 = p; p1[Mt.L.iR] = 1.;
            std::vector<double> mu1 = Mt.Expected(p1);
            const double br = (which >= 3) ? 1e-5 : 0.;
            if (which == 0) { nobs = mu0; I.CB = cbn; label = "smooth parent, background-only Asimov set (CB nominal)"; }
            else {
                std::mt19937_64 rng(0x7c15ULL + 101ULL * which);
                nobs.resize(mu0.size());
                for (size_t b = 0; b < mu0.size(); ++b) {
                    const double m = mu0[b] + br / x17m::kBRdir * (mu1[b] - mu0[b]);
                    std::poisson_distribution<long long> P(std::max(m, 0.));
                    nobs[b] = m > 0 ? (double)P(rng) : 0.;
                }
                I.CB = cbr;
                label = Form("smooth-parent toy at BR = %.0e (CB replica %d)", br, which);
            }
        } else {
            // a simulation-study toy at BR = 0
            TFile* fk = TFile::Open(x17m::kToyKey, "READ");
            TTree* tk = (TTree*)fk->Get("key");
            int id; double br;
            tk->SetBranchAddress("id", &id); tk->SetBranchAddress("br", &br);
            int found = -1, cnt = 0;
            for (Long64_t i = 0; i < tk->GetEntries(); ++i) { tk->GetEntry(i); if (br == 0) { if (cnt == which - 5) { found = id; break; } ++cnt; } }
            fk->Close(); delete fk;
            TFile* fx = TFile::Open(x17m::kToyFile, "READ");
            TH2* h = (TH2*)fx->Get(Form("toy_%d_m4mee", found));
            nobs = x17m::ToyRebin(h, to);
            fx->Close(); delete fx;
            std::mt19937_64 rng(0x5eedULL + found);
            for (size_t b = 0; b < nobs.size(); ++b) { std::poisson_distribution<long long> P(cbn[b]); nobs[b] += cbn[b] > 0 ? (double)P(rng) : 0.; }
            I.CB = cbr;
            label = Form("simulation-study toy %d at BR = 0 (+ CB-like)", found);
        }
    } else {
        // ---------------------------------------------------------------- data-scale projection (x17ab_model.C)
        TFile* f = TFile::Open(alldefs::kHistFile, "READ");
        x17m::HistOpt o;
        if (!x17m::BuildFromHists(f, I, o, alldefs::kNbToRaw)) return;
        f->Close(); delete f;
        TFile* fm = TFile::Open(Form("%s/x17ab_model.root", kOut), "READ");
        TH2* ha = fm ? (TH2*)fm->Get("asimov_proj") : nullptr;
        if (!ha) { printf("x17tc: run x17ab_model.C first\n"); return; }
        nobs.assign(I.G.n(), 0.);
        for (int j = 0; j < I.G.ny; ++j) for (int i = 0; i < I.G.nx; ++i) nobs[I.G.idx(i, j)] = ha->GetBinContent(i + 1, j + 1);
        fm->Close(); delete fm;
        C.a0 = 0.40; C.cTailSig = 0.194;
        label = "data-scale X2D, pre-data projection Asimov set (a_eta 0.40)";
    }
    x17m::Model M(I, C);
    const x17m::Glob g = x17m::GlobOf(C);
    x17m::Limit La = x17m::LimitAsym(M, nobs, g, 1e-5 / x17m::kBRdir, true);
    printf("x17tc: [%d] %s: rhat %.4e +- %.4e; asymptotic UL %.4e (BR %.3e); expected median %.4e (BR %.3e)\n", which, label.c_str(), La.rhat, La.rhatErr, La.ul,
           La.ul * x17m::kBRdir, La.exp[2], La.exp[2] * x17m::kBRdir);
    const double fr[8] = {0.6, 0.75, 0.85, 0.95, 1.05, 1.15, 1.3, 1.6};
    std::vector<x17m::ToyCLsPoint> pts;
    FILE* tf = fopen(Form("%s/toycls_%d.txt", kOut, which), "w");
    fprintf(tf, "# x17ab_toycls.C data set %d: %s; %d toys per point and hypothesis\n", which, label.c_str(), ntoys);
    fprintf(tf, "# mu(r)       BR        q_obs   CLs+b    CLb     CLs    +-   | asymptotic CLs\n");
    for (int k = 0; k < 8; ++k) {
        const double mu = fr[k] * La.ul;
        x17m::ToyCLsPoint T = x17m::ToyCLs(M, nobs, g, La.prof, mu, ntoys, 0xC15ULL * (which + 1) + 7919ULL * k);
        // asymptotic CLs at the same mu
        x17m::Profile PA; PA.zero = M.Fit(La.asimov, 0., false, g, &La.prof.zero.p); PA.hat = PA.zero; PA.rhat = 0; PA.Fref = PA.zero.F;
        const double qA = x17m::QTilde(M, La.asimov, g, PA, mu);
        const double clsA = x17m::CLsAsym(T.qobs, qA);
        pts.push_back(T);
        fprintf(tf, "%.5e %.4e %8.4f %7.4f %7.4f %7.4f %6.4f | %7.4f\n", mu, mu * x17m::kBRdir, T.qobs, T.clsb, T.clb, T.cls, T.clsErr, clsA);
        printf("x17tc: [%d] mu %.4e: q %.3f CLs+b %.4f CLb %.4f CLs %.4f +- %.4f (asymptotic %.4f), %.0f s\n", which, mu, T.qobs, T.clsb, T.clb, T.cls, T.clsErr, clsA,
               sw.RealTime());
        sw.Continue();
    }
    const double ulT = x17m::ToyUL(pts);
    fprintf(tf, "SUMMARY %d asymptotic %.5e toys %.5e ratio %.4f  (BR: %.4e / %.4e)\n", which, La.ul, ulT, ulT / La.ul, La.ul * x17m::kBRdir, ulT * x17m::kBRdir);
    fclose(tf);
    printf("x17tc: [%d] UL asymptotic %.4e, toys %.4e, ratio %.3f\n", which, La.ul, ulT, ulT / La.ul);
}
