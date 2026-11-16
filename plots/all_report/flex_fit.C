// ============================================================================
// flex_fit.C  --  package flex of report_all: flexible SMASH description of the data, fits and variants
// ============================================================================
// What it does (plain words):
//   1. Checks: the fit bins built from flex_data.root / flex_sim.root equal hists_v1.root F10 E2 (data N, CB, templates), and
//      the Part III F-main configuration of this library reproduces wp5 results.txt (a_eta 0.382284, ...).
//   2. "SMASH as is": every component at its SMASH cross section (all scales 1, k_E2 = 1): chi2 only.
//   3. The flexible fit: N = k CB + sum_i s_i T_i with one free scale per component, on F10 E2 (5 M(pi+pi-) slices,
//      M(ee) >= 10 MeV, M4 0.40-0.70). Ladder of models (no tilt / common tilt / per-slice tilts; 7 components / omega-rho0
//      split of G_D and G_C); the nominal one is chosen with the D11 rule of Part III (add a term only if the hybrid GLS chi2,
//      which includes the template MC noise, drops by more than 4).
//   4. Variants of the chosen model: M(ee) >= 10 MeV split into 10-20 and >= 20 MeV; E1 added with a free CB-like term k_E1
//      or with a free fake-lepton term (T3); M4 range to 1.00; ranges 0.42-0.70 and 0.45-0.70; raw background templates.
//   Every fit: value, stat (sandwich: event groups + CB replicas + constraints), MC (100 simulation replicas), correlations,
//   Pearson / Baker-Cousins / hybrid GLS chi2, window contributions.
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/flex/flex_results.txt (machine-readable), flex_fit.root (per-bin parts for
//         the figures), flex_fit.log; docs/latex/tables_all/flex_*.txt
// Run (from /home/przygoda/HADES/FAT/FAT; about 3 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_flex", kTRUE)' 'plots/all_report/flex_fit.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "flex_lib.h"

namespace flexfit {
using namespace flex;

FILE* gRes = nullptr;
FILE* gLog = nullptr;
void L(const char* fmt, ...) __attribute__((format(printf, 1, 2)));
void L(const char* fmt, ...) {
    va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
    va_start(a, fmt); vfprintf(gLog, fmt, a); va_end(a);
}

struct Config {
    std::string tag, text;
    Layout lay;
    ModelOpt o;
    bool p3 = false, rawBg = false, asis = false;
};
struct Out {
    Config cfg;
    Model M;
    Bins B;
    Templ T;
    BinT BT;
    FitRes F;
    RepFits R;
    Summary S;
    GLSRes G;
    std::vector<SmoothInfo> sinfo;
    bool ok = false;
};

const std::vector<int> kK7 = {cC2, cC3, cC4, cC5, cGeta, cGD, cGC};
const std::vector<int> kK9 = {cC2, cC3, cC4, cC5, cGeta, cGDom, cGDnom, cGCom, cGCnom};
const std::vector<int> kK8 = {cC2, cC3, cC4, cC5, cGeta, cGDom, cGDnom, cGC};

Layout LayE2(double lo = 0.40, double hi = 0.70) { Layout L; L.fc = {ClsE2()}; L.mlo = lo; L.mhi = hi; return L; }

std::map<std::string, Bins> gBinsCache;
const Bins& GetBins(const Layout& L) {
    std::string key = Form("%.4f_%.4f", L.mlo, L.mhi);
    for (const auto& f : L.fc) key += "_" + f.tag;
    if (!gBinsCache.count(key)) { Bins B; MakeBins(L, B); gBinsCache[key] = B; }
    return gBinsCache[key];
}
// BB-lite / raw switch for the background (raw background templates with their variance)
BinT BinTemplRaw(const Bins& B, const Templ& T, const Model& M, bool rawBg) {
    BinT o = BinTempl(B, T, M);
    if (!rawBg) return o;
    const std::vector<CompDef> D = AllComps();
    for (int c : M.Comps())
        if (D[c].smooth)
            for (int i = 0; i < B.n; ++i) { o.t[c][i] = T.C(c, B.h[i])[B.x[i]]; o.v[c][i] = T.V(c, B.h[i])[B.x[i]]; }
    return o;
}

Out Run(SimSrc& S, const Config& c, bool reps = true) {
    TStopwatch sw;
    Out o; o.cfg = c;
    bool hasE1 = false, hasE2 = false;
    for (const auto& f : c.lay.fc) { if (f.kidx == 0) hasE1 = true; else hasE2 = true; }
    o.M = c.p3 ? MakePartIII(hasE1, hasE2) : MakeFlex(c.o, hasE1, hasE2, c.tag, c.text);
    o.M.tag = c.tag; o.M.text = c.text;
    o.B = GetBins(c.lay);
    o.T = MakeTempl(S, c.lay, o.M.Comps(), -1, &o.sinfo);
    o.BT = BinTemplRaw(o.B, o.T, o.M, c.rawBg);
    if (c.asis) {
        // all scales at SMASH (1), k = 1: no fit
        FitRes F;
        F.th = o.M.start; F.n = o.B.n; F.nfree = 0;
        Eval(o.B, o.BT, o.M, F.th, F.mu, nullptr, &F.V, &F.Vd);
        for (int i = 0; i < o.B.n; ++i) {
            const double m = std::max(F.mu[i], 1e-9), y = o.B.y[i];
            F.chi2BC += 2. * (m - y + (y > 0 ? y * std::log(y / m) : 0.));
            F.chi2P += (y - F.mu[i]) * (y - F.mu[i]) / F.V[i];
        }
        F.ndf = o.B.n;
        F.Cstat = F.Cn = F.Ccb = F.C = Zero(o.M.np(), o.M.np());
        F.ok = true;
        o.F = F;
        if (reps) { RepSet RS = MakeReplicaSet(S, c.lay, o.B, o.M); o.G = HybridGLS(o.B, o.BT, o.M, o.F, RS); }
        o.S = Summarise(o.M, o.F, nullptr);
        o.ok = true;
        L("[%s] as is: n %d Pearson %.1f BC %.1f hybridGLS %.1f (data-only GLS %.1f) (%.1f s)\n", c.tag.c_str(), o.B.n, F.chi2P, F.chi2BC,
          o.G.chi2, o.G.chi2data, sw.RealTime());
        return o;
    }
    o.F = Fit(o.B, o.BT, o.M, nullptr, true, &o.B.eg);
    if (!o.F.ok) { L("[%s] FIT FAILED (nit %d)\n", c.tag.c_str(), o.F.nit); return o; }
    if (reps) {
        RepSet RS = MakeReplicaSet(S, c.lay, o.B, o.M);
        if (c.rawBg) for (size_t r = 0; r < RS.bt.size(); ++r) RS.bt[r] = BinTemplRaw(o.B, RS.T[r], o.M, true);
        o.R = RunReplicas(o.B, o.M, RS, o.F.th);
        o.G = HybridGLS(o.B, o.BT, o.M, o.F, RS, !c.rawBg);
        o.S = Summarise(o.M, o.F, &o.R);
    } else o.S = Summarise(o.M, o.F, nullptr);
    o.ok = true;
    L("[%s] %s\n  n %d free %d ndf %d nit %d | Pearson %.1f (%.2f) BC %.1f | hybrid GLS %.1f (%.2f), data-only GLS %.1f | reps ok %d (%.1f s)\n",
      c.tag.c_str(), c.text.c_str(), o.B.n, o.F.nfree, o.F.ndf, o.F.nit, o.F.chi2P, o.F.chi2P / o.F.ndf, o.F.chi2BC, o.G.chi2,
      o.G.chi2 / std::max(1, o.F.ndf), o.G.chi2data, o.R.nok, sw.RealTime());
    for (int k = 0; k < o.M.np(); ++k)
        L("    %-10s %9.4f  stat %.4f  MC %.4f  tot %.4f  (rel %.3f)  shift %+.4f%s\n", o.M.names[k].c_str(), o.S.v[k], o.S.eStat[k], o.S.eMC[k],
          o.S.eTot[k], o.S.v[k] != 0 ? o.S.eTot[k] / std::fabs(o.S.v[k]) : 0., o.S.shift[k], o.M.fixed[k] ? "  [fixed]" : "");
    return o;
}

// window sums (bin centre inside the window) of the model parts, raw rows
const int kNWs = 5;
const char* const kWsTag[kNWs] = {"B", "P", "CRU", "Fit", "U70"};
const double kWsLo[kNWs] = {0.4550, 0.53364, 0.565, 0.40, 0.70}, kWsHi[kNWs] = {0.5225, 0.55416, 0.700, 0.70, 1.00};
void WriteOut(const Out& o) {
    const Model& M = o.M;
    const int np = M.np();
    fprintf(gRes, "\nfit %s ok %d\ntext %s %s\n", o.cfg.tag.c_str(), o.ok ? 1 : 0, o.cfg.tag.c_str(), o.cfg.text.c_str());
    if (!o.ok) return;
    std::string fcs;
    for (const auto& f : o.cfg.lay.fc) fcs += f.tag + " ";
    fprintf(gRes, "layout %s classes %s range %.4f %.4f nbins %d\n", o.cfg.tag.c_str(), fcs.c_str(), o.cfg.lay.mlo, o.cfg.lay.mhi, o.B.n);
    fprintf(gRes, "gof %s nbins %d nfree %d ndf %d chi2P %.3f chi2BC %.3f chi2GLShyb %.3f chi2GLSdata %.3f nrep_ok %d asis %d\n", o.cfg.tag.c_str(), o.B.n,
            o.F.nfree, o.F.ndf, o.F.chi2P, o.F.chi2BC, o.G.chi2, o.G.chi2data, o.R.nok, o.cfg.asis ? 1 : 0);
    auto row = [&](const char* key, const Vec& v) {
        fprintf(gRes, "%s %s", key, o.cfg.tag.c_str());
        for (int k = 0; k < np; ++k) fprintf(gRes, " %.8g", v[k]);
        fprintf(gRes, "\n");
    };
    fprintf(gRes, "names %s", o.cfg.tag.c_str()); for (auto& n : M.names) fprintf(gRes, " %s", n.c_str()); fprintf(gRes, "\n");
    fprintf(gRes, "fixed %s", o.cfg.tag.c_str()); for (int k = 0; k < np; ++k) fprintf(gRes, " %d", M.fixed[k]); fprintf(gRes, "\n");
    row("values", o.S.v); row("err_N", o.S.eN); row("err_CB", o.S.eCB); row("err_stat", o.S.eStat); row("err_mc", o.S.eMC); row("err_total", o.S.eTot);
    row("mc_shift", o.S.shift);
    fprintf(gRes, "cor_total %s %d\n", o.cfg.tag.c_str(), np);
    for (int k = 0; k < np; ++k) { for (int l = 0; l < np; ++l) fprintf(gRes, " %7.4f", o.S.corTot[k][l]); fprintf(gRes, "\n"); }
    fprintf(gRes, "cov_total %s %d\n", o.cfg.tag.c_str(), np);
    for (int k = 0; k < np; ++k) { for (int l = 0; l < np; ++l) fprintf(gRes, " %.6g", o.S.covTot[k][l]); fprintf(gRes, "\n"); }
    fprintf(gRes, "cov_stat %s %d\n", o.cfg.tag.c_str(), np);
    for (int k = 0; k < np; ++k) { for (int l = 0; l < np; ++l) fprintf(gRes, " %.6g", o.S.covStat[k][l]); fprintf(gRes, "\n"); }
    // smoothing information
    const std::vector<CompDef> D = AllComps();
    for (int c : M.Comps()) {
        if (!D[c].smooth) continue;
        for (int h = 0; h < o.T.nh; ++h) {
            const SmoothInfo& s = o.sinfo[c * o.T.nh + h];
            fprintf(gRes, "smooth %s %s %s neff %.1f order %d scaled %d chi2 %.1f nb %d", o.cfg.tag.c_str(), D[c].name.c_str(), o.cfg.lay.HTag(h).c_str(),
                    s.neff[0], s.order[0], s.scaled[0] ? 1 : 0, s.chi2[0], s.nb[0]);
            if (s.nb[1] > 0) fprintf(gRes, " | hi neff %.1f order %d scaled %d chi2 %.1f nb %d", s.neff[1], s.order[1], s.scaled[1] ? 1 : 0, s.chi2[1], s.nb[1]);
            fprintf(gRes, "\n");
        }
    }
    // window contributions at the fitted scales and at SMASH (raw rows), per class k index
    Vec mu; std::vector<Vec> parts, parts1;
    Eval(o.B, o.BT, M, o.F.th, mu, nullptr, nullptr, nullptr, &parts);
    Eval(o.B, o.BT, M, M.start, mu, nullptr, nullptr, nullptr, &parts1);
    for (int w = 0; w < kNWs; ++w) {
        double N = 0, cb = 0, tot = 0;
        std::map<int, double> s, s1;
        for (int i = 0; i < o.B.n; ++i) {
            const double xc = 0.5 * (o.B.lo[i] + o.B.hi[i]);
            if (xc < kWsLo[w] || xc >= kWsHi[w]) continue;
            N += o.B.y[i]; cb += parts[kNCompF][i];
            for (int c : M.Comps()) { s[c] += parts[c][i]; s1[c] += parts1[c][i]; tot += parts[c][i]; }
        }
        if (N == 0) continue;
        fprintf(gRes, "win %s %s N %.1f kCB %.1f model %.1f", o.cfg.tag.c_str(), kWsTag[w], N, cb, cb + tot);
        for (int c : M.Comps()) fprintf(gRes, " %s %.1f/%.1f", D[c].name.c_str(), s[c], s1[c]);
        fprintf(gRes, "\n");
    }
    fflush(gRes);
}

// histograms of the parts for the figures (one TH1D per fit histogram and part, X axis bins of the layout)
void WriteHists(TFile* f, const Out& o, const char* dirName, const Vec& th) {
    if (!o.ok || th.empty()) return;
    TDirectory* d = f->mkdir(dirName);
    const std::vector<double> E = XEdges();
    Vec mu, V; std::vector<Vec> parts;
    Eval(o.B, o.BT, o.M, th, mu, nullptr, &V, nullptr, &parts);
    const std::vector<CompDef> D = AllComps();
    for (int h = 0; h < o.cfg.lay.nh(); ++h) {
        auto mk = [&](const std::string& n) { TH1D* x = new TH1D((n + "_" + o.cfg.lay.HTag(h)).c_str(), "", kNX, E.data()); x->SetDirectory(nullptr); x->Sumw2(); return x; };
        TH1D* hN = mk("N"); TH1D* hCB = mk("kCB"); TH1D* hMu = mk("mu"); TH1D* hEG = mk("vEG");
        std::map<int, TH1D*> hp;
        for (int c : o.M.Comps()) hp[c] = mk(D[c].name);
        for (int i = 0; i < o.B.n; ++i) {
            if (o.B.h[i] != h) continue;
            const int x = o.B.x[i];
            hN->SetBinContent(x, o.B.y[i]); hN->SetBinError(x, std::sqrt(std::max(0., o.B.eg[(size_t)i * o.B.n + i])));
            hCB->SetBinContent(x, parts[kNCompF][i]); hCB->SetBinError(x, std::sqrt(o.B.vcb[i]));
            hMu->SetBinContent(x, mu[i]); hMu->SetBinError(x, std::sqrt(V[i]));
            hEG->SetBinContent(x, o.B.eg[(size_t)i * o.B.n + i]);
            for (int c : o.M.Comps()) hp[c]->SetBinContent(x, parts[c][i]);
        }
        d->WriteTObject(hN); d->WriteTObject(hCB); d->WriteTObject(hMu); d->WriteTObject(hEG);
        for (auto& p : hp) d->WriteTObject(p.second);
    }
    // parameter values and names
    if (!o.ok || th.empty()) return;
    TH1D* hv = new TH1D("theta", "", o.M.np(), 0, o.M.np());
    hv->SetDirectory(nullptr);
    for (int k = 0; k < o.M.np(); ++k) { hv->GetXaxis()->SetBinLabel(k + 1, o.M.names[k].c_str()); hv->SetBinContent(k + 1, th[k]); }
    d->WriteTObject(hv);
}

// derived quantity: value, stat (delta method with Cstat), MC (replica refits)
struct DQ { double v = 0, eStat = 0, eMC = 0, eTot = 0; };
DQ Derived(const Out& o, const std::function<double(const Vec&)>& f) {
    DQ d;
    d.v = f(o.F.th);
    const int np = o.M.np();
    Vec g(np, 0.);
    for (int k = 0; k < np; ++k) {
        if (o.M.fixed[k]) continue;
        const double h = std::max(1e-7, 1e-4 * std::sqrt(std::max(o.F.Cstat[k][k], 1e-12)));
        Vec a = o.F.th, b = o.F.th; a[k] += h; b[k] -= h;
        g[k] = (f(a) - f(b)) / (2 * h);
    }
    d.eStat = std::sqrt(std::max(0., Quad(o.F.Cstat, g)));
    double s = 0, q = 0; int n = 0;
    for (const auto& r : o.R.fit) if (r.ok) { const double v = f(r.th); s += v; q += v * v; ++n; }
    if (n > 1) { const double m = s / n; d.eMC = std::sqrt(std::max(0., (q - n * m * m) / (n - 1))); }
    d.eTot = std::sqrt(d.eStat * d.eStat + d.eMC * d.eMC);
    return d;
}

}  // namespace flexfit

void flex_fit() {
    using namespace flex;
    using namespace flexfit;
    TStopwatch total;
    gSystem->mkdir(kOut, kTRUE);
    gLog = fopen(Form("%s/flex_fit.log", kOut), "w");
    gRes = fopen(Form("%s/flex_results.txt", kOut), "w");
    fprintf(gRes, "# flex_results.txt -- plots/all_report/flex_fit.C (package flex: flexible SMASH description of the data)\n");
    fprintf(gRes, "# S0, isBest >= 0, L = 6079.2 nb^-1, RAW rows; templates = SMASH sigma x PDG BR x V3 at M4' = 0.99615 M4 (nb x %.4f)\n", kNbToRaw);
    fprintf(gRes, "# model: mu = k_e CB + sum_t s_t T_t [x c_j (1 + b x) for background terms]; scales s = 1 is SMASH as is\n");
    SimSrc S;
    if (!S.Open()) { printf("flex: cannot open %s\n", kSimFile); return; }

    // ------------------------------------------------------------------ 1. checks
    {
        const Layout L0 = LayE2();
        const Bins& B = GetBins(L0);
        TFile* fh = TFile::Open(kHistFile);
        double dy = 0, dcb = 0, dc2 = 0, dgd = 0, dvcb = 0;
        for (int i = 0; i < B.n; ++i) {
            const std::string hn = F10Name(B.j[i], 1);
            TH1* n = allio::Data(fh, "S0", "all", "N", hn);
            TH1* cb = allio::Data(fh, "S0", "all", "CB", hn);
            TH1* c2 = allio::Sim(fh, "S0", "nom", kC2, hn);
            TH1* gd = allio::Sim(fh, "S0", "nom", kGD, hn);
            dy = std::max(dy, std::fabs(n->GetBinContent(B.x[i]) - B.y[i]));
            dcb = std::max(dcb, std::fabs(cb->GetBinContent(B.x[i]) - B.cb[i]));
            delete n; delete cb;
            Vec c, v;
            SumParts(S, {kC2}, B.j[i], ClsE2(), -1, c, v);
            dc2 = std::max(dc2, std::fabs(c[B.x[i]] - c2->GetBinContent(B.x[i]) * kNbToRaw));
            SumParts(S, {kGD}, B.j[i], ClsE2(), -1, c, v);
            dgd = std::max(dgd, std::fabs(c[B.x[i]] - gd->GetBinContent(B.x[i]) * kNbToRaw));
            delete c2; delete gd;
        }
        // CB variance against the wp5 definition (replica covariance of F10 CBrep, 1/n)
        for (int j = 0; j < kNSlice; ++j) {
            TH2* r = (TH2*)fh->Get(Form("data/S0/all/CBrep_%s", F10Name(j, 1).c_str()));
            for (int i = 0; i < B.n; ++i) {
                if (B.j[i] != j) continue;
                double m = 0, q = 0;
                for (int k = 0; k < kNRepC; ++k) { const double x = r->GetBinContent(B.x[i], k + 1); m += x; q += x * x; }
                m /= kNRepC; q = q / kNRepC - m * m;
                dvcb = std::max(dvcb, std::fabs(q - B.vcb[i]) / std::max(1e-9, q));
            }
        }
        fh->Close();
        L("check: F10 E2 0.40-0.70 (%d bins) against hists_v1.root: max|dN| %.3g, max|dCB| %.3g, max rel d var(CB) %.3g, max|dC2| %.3g, max|dGD| %.3g raw rows\n",
          B.n, dy, dcb, dvcb, dc2, dgd);
        fprintf(gRes, "check hists_v1 F10E2 bins %d dN %.3g dCB %.3g dvarCB_rel %.3g dC2 %.3g dGD %.3g\n", B.n, dy, dcb, dvcb, dc2, dgd);
    }
    std::vector<Out> outs;
    TFile* fr = TFile::Open(Form("%s/flex_fit.root", kOut), "RECREATE");
    // Part III reproduction
    {
        Config c; c.tag = "P3"; c.text = "Part III F-main reproduced (E2, kappa = 1 +- 0.2, g = 1 +- 0.3): check against wp5 results.txt"; c.lay = LayE2(); c.p3 = true;
        Out o = Run(S, c, true);
        const double ref[12] = {0.38228397, 0.83233498, 1.3152029, 1.7167655, 1.0922773, 1.1046459, 1.090109, 1.1529077, 0.98427253, -0.011082838, 1.9469467, 0.98181167};
        const char* rn[12] = {"a_eta", "kappa", "rho_d", "g", "c_1", "c_2", "c_3", "c_4", "c_5", "b", "kappa_B", "k_E2"};
        double dmax = 0;
        for (int q = 0; q < 12; ++q) { const int k = o.M.Idx(rn[q]); dmax = std::max(dmax, std::fabs(o.F.th[k] - ref[q])); }
        L("check: Part III F-main reproduced: max |theta - wp5| = %.3g (a_eta %.6f vs 0.382284; stat %.5f vs 0.016925; MC %.5f vs 0.039629)\n", dmax,
          o.F.th[o.M.Idx("a_eta")], o.S.eStat[o.M.Idx("a_eta")], o.S.eMC[o.M.Idx("a_eta")]);
        fprintf(gRes, "check P3 max_dtheta %.3g a_eta %.6f stat %.6f mc %.6f\n", dmax, o.F.th[o.M.Idx("a_eta")], o.S.eStat[o.M.Idx("a_eta")], o.S.eMC[o.M.Idx("a_eta")]);
        WriteOut(o);
        WriteHists(fr, o, "P3", o.F.th);
        // Part III-equivalent scales of the flexible components (weighted over the fit bins)
        const int ia = o.M.Idx("a_eta"), ik = o.M.Idx("kappa"), ir = o.M.Idx("rho_d"), ig = o.M.Idx("g"), ikb = o.M.Idx("kappa_B");
        Vec mu; std::vector<Vec> parts, parts1;
        Eval(o.B, o.BT, o.M, o.F.th, mu, nullptr, nullptr, nullptr, &parts);
        Vec th1 = o.F.th;
        for (int j = 0; j < kNSlice; ++j) th1[o.M.iC[j]] = 1.;
        th1[o.M.iB] = 0.; th1[ikb] = 1.;
        Eval(o.B, o.BT, o.M, th1, mu, nullptr, nullptr, nullptr, &parts1);
        double sgd = 0, sgd1 = 0, sgc = 0, sgc1 = 0;
        for (int i = 0; i < o.B.n; ++i) { sgd += parts[cGD][i]; sgd1 += parts1[cGD][i]; sgc += parts[cGC][i]; sgc1 += parts1[cGC][i]; }
        fprintf(gRes, "p3equiv C2 %.6f C3 %.6f C4 %.6f C5 %.6f Geta %.6f GD %.6f GC %.6f k_E2 %.6f\n", o.F.th[ia], o.F.th[ia] * o.F.th[ik], o.F.th[ia] * o.F.th[ik],
                o.F.th[ia] * o.F.th[ir], o.F.th[ia] * o.F.th[ig], sgd / sgd1, sgc / sgc1, o.F.th[o.M.Idx("k_E2")]);
        // errors of the equivalents (stat, MC) by the delta method / replicas
        auto eq = [&](int which) {
            return [=, &o](const Vec& th) {
                switch (which) {
                    case 0: return th[ia];
                    case 1: return th[ia] * th[ik];
                    case 2: return th[ia] * th[ir];
                    case 3: return th[ia] * th[ig];
                }
                return 0.;
            };
        };
        const char* eqn[4] = {"C2", "C3C4", "C5", "Geta"};
        for (int q = 0; q < 4; ++q) { DQ d = Derived(o, eq(q)); fprintf(gRes, "p3equiv_err %s %.6f stat %.6f mc %.6f tot %.6f\n", eqn[q], d.v, d.eStat, d.eMC, d.eTot); }
        // GD and GC equivalents with errors
        auto bgEq = [&](int comp) {
            return [=, &o](const Vec& th) {
                Vec m1; std::vector<Vec> p, p1; Vec t1 = th;
                for (int j = 0; j < kNSlice; ++j) t1[o.M.iC[j]] = 1.;
                t1[o.M.iB] = 0.; t1[ikb] = 1.;
                Eval(o.B, o.BT, o.M, th, m1, nullptr, nullptr, nullptr, &p);
                Eval(o.B, o.BT, o.M, t1, m1, nullptr, nullptr, nullptr, &p1);
                double a = 0, b = 0;
                for (int i = 0; i < o.B.n; ++i) { a += p[comp][i]; b += p1[comp][i]; }
                return a / b;
            };
        };
        for (int comp : {cGD, cGC}) { DQ d = Derived(o, bgEq(comp)); fprintf(gRes, "p3equiv_err %s %.6f stat %.6f mc %.6f tot %.6f\n", comp == cGD ? "GD" : "GC", d.v, d.eStat, d.eMC, d.eTot); }
        outs.push_back(o);
    }

    // ------------------------------------------------------------------ 2. SMASH as is
    {
        Config c; c.tag = "ASIS"; c.text = "SMASH as is: every scale 1, k_E2 = 1 (no fit)"; c.lay = LayE2(); c.asis = true; c.o.comps = kK7;
        Out o = Run(S, c, true);
        WriteOut(o);
        WriteHists(fr, o, "ASIS", o.M.start);
        outs.push_back(o);
        Config c9 = c; c9.tag = "ASIS9"; c9.o.comps = kK9;
        Out o9 = Run(S, c9, true);
        WriteOut(o9);
        WriteHists(fr, o9, "ASIS9", o9.M.start);
    }

    // ------------------------------------------------------------------ 3. ladder on F10 E2 0.40-0.70
    struct Lad { const char* tag; const char* text; std::vector<int> comps; int tilt; bool cj; };
    const std::vector<Lad> lad = {
        {"F7", "7 free scales (C2, C3, C4, C5, G_eta, G_D, G_C), no tilt", kK7, 0, false},
        {"F7b", "7 free scales + common tilt (1 + b x) on G_D, G_C", kK7, 1, false},
        {"F7bj", "7 free scales + per-slice tilts (1 + b_j x) on G_D, G_C", kK7, 2, false},
        {"F8", "8 free scales (G_D split omega/rho0 vs other), no tilt", kK8, 0, false},
        {"F8b", "8 free scales (G_D split) + common tilt", kK8, 1, false},
        {"F8bj", "8 free scales (G_D split) + per-slice tilts", kK8, 2, false},
        {"F9", "9 free scales (G_D and G_C split omega/rho0 vs other), no tilt", kK9, 0, false},
        {"F9b", "9 free scales (G_D and G_C split) + common tilt", kK9, 1, false},
        {"F9bj", "9 free scales (G_D and G_C split) + per-slice tilts", kK9, 2, false},
        {"F7cj", "Part III background c_j (1 + b x)[G_D + s_GC G_C] (s_GD = 1) with free C2, C3, C4, C5, G_eta", kK7, 1, true},
    };
    std::map<std::string, Out> ladOut;
    for (const auto& l : lad) {
        Config c; c.tag = l.tag; c.text = l.text; c.lay = LayE2(); c.o.comps = l.comps; c.o.tilt = l.tilt; c.o.perSliceNorm = l.cj;
        Out o = Run(S, c, true);
        WriteOut(o);
        WriteHists(fr, o, l.tag, o.F.th);
        ladOut[l.tag] = o;
    }
    // ------------------------------------------------------------------ model choice
    // D11 rule of Part III applied to the flexible model (hybrid GLS chi2 at the fit, which includes the template MC noise):
    // a term group is added only if chi2 drops by more than 4. Order: the physics splits first (1) omega/rho0 split of G_D,
    // (2) split of G_C, then the shape freedom (3) common tilt, (4) per-slice tilts. The tilt-first order gives the
    // ALTERNATIVE model. The plain rule ends at F9bj, where s_GC_om = -5 +- 2.7 (unphysical) and s_GD_om is undetermined:
    // splits and tilts are degenerate. Admissibility conditions are therefore added to the D11 rule (documented in the note):
    //   (b) every new scale is determined to better than 50 % (sigma_tot / |s| < 0.5);
    //   (c) no scale that was determined to better than 50 % becomes worse than 50 % (no new degeneracy; s_GD -> s_GD_nom);
    //   (d) every scale stays physical: s > -2 sigma_tot.
    auto chi = [&](const std::string& t) { return (ladOut.count(t) && ladOut[t].ok) ? ladOut[t].G.chi2 : 1e30; };
    fprintf(gRes, "\n# D11 rule (hybrid GLS chi2 at the fit, add a term group only if dchi2 > 4) + admissibility (b, c, d); grid:\n");
    for (const auto& l : lad) fprintf(gRes, "grid %s ok %d chi2hyb %.2f chi2P %.2f nfree %d\n", l.tag, ladOut[l.tag].ok ? 1 : 0, chi(l.tag), ladOut[l.tag].F.chi2P, ladOut[l.tag].F.nfree);
    auto relErr = [](const Out& o, int k) { return o.S.v[k] != 0 ? o.S.eTot[k] / std::fabs(o.S.v[k]) : 1e9; };
    auto step = [&](std::string& cur, const std::string& cand, const char* what) {
        const double d = chi(cur) - chi(cand);
        std::string why;
        bool add = d > 4 && chi(cand) < 1e29;
        if (!(d > 4)) why += " (a) dchi2 <= 4;";
        if (chi(cand) >= 1e29) why += " fit failed;";
        if (chi(cand) < 1e29 && chi(cur) < 1e29) {
            const Out& a = ladOut[cur];
            const Out& b = ladOut[cand];
            for (int k = 0; k < b.M.np(); ++k) {
                if (b.M.fixed[k] || b.M.names[k].rfind("s_", 0) != 0) continue;
                const bool isNew = a.M.Idx(b.M.names[k]) < 0 && !(b.M.names[k].size() > 4 && b.M.names[k].substr(b.M.names[k].size() - 4) == "_nom");
                if (isNew && relErr(b, k) >= 0.5) { add = false; why += Form(" (b) new %s undetermined (%.2f);", b.M.names[k].c_str(), relErr(b, k)); }
                if (b.S.v[k] < -2 * b.S.eTot[k]) { add = false; why += Form(" (d) %s = %.2f +- %.2f unphysical;", b.M.names[k].c_str(), b.S.v[k], b.S.eTot[k]); }
            }
            for (int k = 0; k < a.M.np(); ++k) {
                if (a.M.fixed[k] || a.M.names[k].rfind("s_", 0) != 0 || relErr(a, k) >= 0.5) continue;
                int kb = b.M.Idx(a.M.names[k]);
                if (kb < 0) kb = b.M.Idx(a.M.names[k] + "_nom");
                if (kb >= 0 && relErr(b, kb) >= 0.5) { add = false; why += Form(" (c) %s becomes undetermined (%.2f -> %.2f);", a.M.names[k].c_str(), relErr(a, k), relErr(b, kb)); }
            }
        }
        fprintf(gRes, "d11 %s -> %s (%s) dchi2 %.2f %s%s\n", cur.c_str(), cand.c_str(), what, d, add ? "ADD" : "no:", add ? "" : why.c_str());
        L("D11: %s -> %s (%s) dchi2_hyb %.2f %s%s\n", cur.c_str(), cand.c_str(), what, d, add ? "ADD" : "no:", add ? "" : why.c_str());
        if (add) cur = cand;
    };
    std::string nominal, alternative;
    {
        std::string cur = "F7";
        step(cur, "F8", "omega/rho0 split of G_D");
        if (cur == "F8") step(cur, "F9", "split of G_C");
        const std::string comp = cur;
        step(cur, comp + "b", "common tilt");
        step(cur, comp + "bj", "per-slice tilts");
        nominal = cur;
        std::string c2 = "F7";
        step(c2, "F7b", "[tilt-first] common tilt");
        step(c2, "F7bj", "[tilt-first] per-slice tilts");
        const std::string sfx = c2.substr(2);
        step(c2, "F8" + sfx, "[tilt-first] omega/rho0 split of G_D");
        if (c2 == "F8" + sfx) step(c2, "F9" + sfx, "[tilt-first] split of G_C");
        alternative = c2;
    }
    fprintf(gRes, "nominal %s\nalternative %s\n", nominal.c_str(), alternative.c_str());
    L("NOMINAL flexible model: %s; ALTERNATIVE (tilt-first): %s\n", nominal.c_str(), alternative.c_str());
    for (const std::string& which : {nominal, alternative}) {
        const Out& NO = ladOut[which];
        const Model& M = NO.M;
        auto sc = [&](int c) { return M.scaleOf[c]; };
        struct R { const char* name; int a, b; };
        std::vector<R> rs = {{"s_C3/s_C2 (kappa-like)", sc(cC3), sc(cC2)}, {"s_C4/s_C2 (kappa-like)", sc(cC4), sc(cC2)},
                             {"s_C5/s_C2 (rho_d-like)", sc(cC5), sc(cC2)}, {"s_Geta/s_C2 (g-like)", sc(cGeta), sc(cC2)}};
        if (sc(cGD) >= 0 && sc(cGC) >= 0) rs.push_back({"s_GC/s_GD (kappa_B-like)", sc(cGC), sc(cGD)});
        if (sc(cGDnom) >= 0 && sc(cGC) >= 0) rs.push_back({"s_GC/s_GD_nom (kappa_B-like)", sc(cGC), sc(cGDnom)});
        if (sc(cGDom) >= 0 && sc(cGDnom) >= 0) rs.push_back({"s_GD_om/s_GD_nom", sc(cGDom), sc(cGDnom)});
        for (const auto& r : rs) {
            if (r.a < 0 || r.b < 0) continue;
            const int a = r.a, b = r.b;
            DQ d = Derived(NO, [=](const Vec& th) { return th[a] / th[b]; });
            fprintf(gRes, "ratio %s %s %.6f stat %.6f mc %.6f tot %.6f\n", which.c_str(), r.name, d.v, d.eStat, d.eMC, d.eTot);
            L("ratio %-6s %-28s %.4f +- %.4f (stat) +- %.4f (MC)\n", which.c_str(), r.name, d.v, d.eStat, d.eMC);
        }
        // content-weighted mean scales of sums over the fit bins (the data may fix a sum better than its parts)
        const std::vector<std::pair<std::string, std::vector<int>>> sums = {{"bump_C2+C3", {cC2, cC3}}, {"peak_C4+C5", {cC4, cC5}}, {"conv_C3+C4", {cC3, cC4}},
                                                                            {"eta_all", {cC2, cC3, cC4, cC5, cGeta}}};
        for (const auto& s : sums) {
            std::vector<int> cs = s.second;
            double w1 = 0;
            for (int c : cs) for (int i = 0; i < NO.B.n; ++i) w1 += NO.BT.t[c][i];
            const Bins* pB = &NO.B; const BinT* pT = &NO.BT;
            std::vector<int> idx; for (int c : cs) idx.push_back(M.scaleOf[c]);
            DQ d = Derived(NO, [=](const Vec& th) { double a = 0; for (size_t q = 0; q < cs.size(); ++q) for (int i = 0; i < pB->n; ++i) a += th[idx[q]] * pT->t[cs[q]][i]; return a / w1; });
            fprintf(gRes, "sumscale %s %s %.6f stat %.6f mc %.6f tot %.6f  # content-weighted mean scale over the fit bins\n", which.c_str(), s.first.c_str(), d.v, d.eStat, d.eMC, d.eTot);
            L("sumscale %-6s %-12s %.4f +- %.4f (stat) +- %.4f (MC)\n", which.c_str(), s.first.c_str(), d.v, d.eStat, d.eMC);
        }
    }
    // ------------------------------------------------------------------ 4. variants of the nominal (V_) and alternative (A_) models
    struct Var { std::string tag, text; Layout lay; ModelOpt o; bool rawBg; };
    std::vector<Var> vars;
    for (int q = 0; q < 2; ++q) {
        const std::string mt = q == 0 ? nominal : alternative;
        if (q == 1 && alternative == nominal) break;
        const std::string P = q == 0 ? "V_" : "A_";
        const Lad* nl = nullptr;
        for (const auto& l : lad) if (mt == l.tag) nl = &l;
        auto base = [&]() { ModelOpt o; o.comps = nl->comps; o.tilt = nl->tilt; o.perSliceNorm = nl->cj; return o; };
        const std::string mtx = " [model " + mt + "]";
        {
            Layout l; l.fc = {ClsE2a(), ClsE2b()};
            vars.push_back({P + "E2split", "M(ee) >= 10 MeV split into 10-20 and >= 20 MeV (k_E2 shared)" + mtx, l, base(), false});
        }
        {
            Layout l; l.fc = {ClsE1(), ClsE2()};
            vars.push_back({P + "E1k", "E1 + E2 with a free CB-like term k_E1 in E1" + mtx, l, base(), false});
            ModelOpt o = base(); o.kE1free = false; o.comps.push_back(cT3); o.e1only = {cT3};
            vars.push_back({P + "E1T3", "E1 + E2 with a free fake-lepton term s_T3 T3 in E1 (k_E1 = 1)" + mtx, l, o, false});
            ModelOpt o2 = base(); o2.comps.push_back(cT3); o2.e1only = {cT3};
            vars.push_back({P + "E1kT3", "E1 + E2 with k_E1 and s_T3 (E1) both free" + mtx, l, o2, false});
            Layout l3; l3.fc = {ClsE1(), ClsE2a(), ClsE2b()};
            vars.push_back({P + "E1k3", "E1 + E2a + E2b with a free CB-like term k_E1" + mtx, l3, base(), false});
        }
        vars.push_back({P + "r100", "M4 range 0.40-1.00 (F10 bins + 10 MeV bins to 1.00, smoothing per piece)" + mtx, LayE2(0.40, 1.00), base(), false});
        vars.push_back({P + "r100raw", "M4 range 0.40-1.00 with raw background templates (BB-lite; no smoothing step at 0.70)" + mtx, LayE2(0.40, 1.00), base(), true});
        vars.push_back({P + "r42", "M4 range 0.42-0.70" + mtx, LayE2(0.42, 0.70), base(), false});
        vars.push_back({P + "r45", "M4 range 0.45-0.70" + mtx, LayE2(0.45, 0.70), base(), false});
        vars.push_back({P + "rawBg", "raw (unsmoothed) background templates with BB-lite variance" + mtx, LayE2(), base(), true});
        {
            ModelOpt o = base(); o.comps.push_back(cT3);
            vars.push_back({P + "T3all", "fake-lepton template T3 added with a free scale in E2 (D6 excludes it)" + mtx, LayE2(), o, false});
        }
        if (q == 0) {
            ModelOpt o = base();
            ModelOpt o7 = o; o7.comps = kK7;
            if (o.comps != kK7) vars.push_back({"V_K7r100", "7 components, no tilt, range 0.40-1.00", LayE2(0.40, 1.00), o7, false});
            ModelOpt o9 = o; o9.comps = kK9;
            if (o.comps != kK9) vars.push_back({"V_K9r100", "9 components (omega/rho0 split of G_D and G_C), range 0.40-1.00" + mtx, LayE2(0.40, 1.00), o9, false});
            ModelOpt ob = o; ob.tilt = 2;
            if (o.tilt != 2) vars.push_back({"V_bjr100", "nominal components + per-slice tilts, range 0.40-1.00", LayE2(0.40, 1.00), ob, false});
        }
    }
    for (const auto& v : vars) {
        Config c; c.tag = v.tag; c.text = v.text; c.lay = v.lay; c.o = v.o; c.rawBg = v.rawBg;
        Out o = Run(S, c, true);
        WriteOut(o);
        WriteHists(fr, o, v.tag.c_str(), o.F.th);
        if (v.tag == "V_r100" || v.tag == "V_r100raw" || v.tag == "V_E1k" || v.tag == "V_E1T3" || v.tag == "V_E1k3") {
            Config ca = c; ca.tag = v.tag + "_asis"; ca.asis = true;
            Out oa = Run(S, ca, false);
            WriteHists(fr, oa, ca.tag.c_str(), oa.M.start);
        }
    }
    fr->Close();
    fclose(gRes);
    L("flex_fit done (%.1f s)\n", total.RealTime());
    fclose(gLog);
}
