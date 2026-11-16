// ============================================================================
// wp7_syst.C  --  WP7 systematics (SPEC report_all Sect. 5), prefix wp7: refit of every variant
// ============================================================================
// Written by WP7, 2026-11-14. Library: wp7_lib.h (generalised WP4/WP5/WP6 code).
// For every configuration (nominal and the Sect. 5 variants S1-S9, S11; consistency tests: folds, periods, isBest == 1,
// cor_ok fallback) the full chain is rerun from the rows (data) and hists_v1.root (simulation):
//   F-main (F10 E2, kappa = 1 +- 0.2, g = 1 +- 0.3, k_E2 = 1 +- 0.05; WP5 nominal after the E1 fall-back),
//   F-gold (F5g; kappa_B and a kappa (C3 + C4) from F-main),
//   stage A (WP3 frozen configuration M1, E2, CR-U 0.575-0.700, eta terms fixed to F-main; WP3 smoothing),
//   B4 (E2, extrapolated to all M(ee)) and the peak side bands (all M(ee)) = the MAIN estimators (WP4/WP6),
// with the statistical error of every quantity (fits: sandwich; side-band estimators: N + CB + stage A + F-main).
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/wp7/syst_raw.txt (one line per configuration and quantity),
//         wp7_syst.log, periods.txt (stability test), checks.
// Run (from /home/przygoda/HADES/FAT/FAT; about 5-10 min, < 6 GB):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp7", kTRUE)' 'plots/all_report/wp7_syst.C+O'
// ============================================================================
#include <fstream>
#include <sstream>
#include "TStopwatch.h"
#include "TSystem.h"
#include "wp7_lib.h"

namespace wp7 {
FILE* gLog = nullptr;
FILE* gRaw = nullptr;
FILE* gRepF = nullptr;   // per-replica values of every quantity (syst_rep.txt)
#define P(...) do { printf(__VA_ARGS__); if (gLog) fprintf(gLog, __VA_ARGS__); fflush(stdout); } while (0)

struct QV { double v = NAN, stat = NAN, eN = NAN, eCB = NAN, eA = NAN, eF = NAN, mc = NAN; int nrep = 0; Vec rep; };
struct Res {
    std::string tag, fam, text;
    bool ok = false;
    std::map<std::string, QV> q;
    double chiF = 0, chiA = 0, chiG = 0; int ndfF = 0, ndfA = 0, ndfG = 0;
    Vec thF, thA, thG;
    Mat CF;
    double kBF = 0;
};
struct Cfg {
    std::string tag, fam, text;
    Layout L = StdLayout();
    SelCfg S;
    CBOpt C;
    int src = 0;                     // 0 F10 of sim/<simSel>__<simVar>; 1 G3 of sim/<simSel>__nom
    std::string simSel = "S0", simVar = "nom";
    double simScale = 1.;
    bool fact = false;
    int ordF = 3, ordA = 3;
    bool rawF = false, rawA = false;
    EtaOpt eo;
    BkgOpt boF, boA;
    int useE = 2;
    double fLo = 0.40, fHi = 0.70;
    double aLo = 0.575, aHi = 0.700;
    bool crk = false;
    int wU = kWU575;
    bool aKBfromKappa = false;       // stage A: kappa_B fixed to the F-main kappa (S3 kappa_B = kappa)
    bool extF = false;               // F-main not refitted: (a, kappa, rho, g) and kappa_B given
    Vec thFext; double kBext = 0; Mat CFext;
    bool gold = true;
    bool bkgFixed = false;           // period test: F-main background fixed to the nominal shape, one scale
    const Res* nomRef = nullptr;     // nominal result (for bkgFixed)
    bool mc = false;                 // run the 100-replica chain (S0__nom templates only)
    bool fixFrep = false;            // F-main fixed to the nominal: replica r uses the nominal replica-r F-main
};

std::map<std::string, DataL> gData;
std::map<std::string, Tmpl> gTmpl;
Rows gRows;
TFile* gFH = nullptr;
std::vector<Tmpl> gRep;          // the 100 replicas of S0__nom
std::vector<Tmpl> gRepPre;       // the 100 replicas of pre__nom (same event-group weights as gRep)
std::vector<Vec> gNomRepThF;     // nominal F-main (a, kappa, rho, g) of every replica (NaN if the fit failed)

DataL& GetData(const Cfg& c) {
    const std::string key = c.L.tag + "|" + c.S.Tag() + "|" + c.C.Tag();
    auto it = gData.find(key);
    if (it != gData.end()) return it->second;
    TStopwatch sw;
    DataL& D = gData[key];
    D.L = c.L; D.S = c.S; D.C = c.C;
    BuildData(gRows, D, true);
    P("  [data %s: e+e- rows %ld, like-sign %ld, items %d, %.1f s]\n", key.c_str(), D.nUS, D.nLS, D.nI, sw.RealTime());
    return D;
}
Tmpl& GetTmpl(const Cfg& c) {
    const std::string key = Form("%d|%s|%s|%s|%.6f|%d", c.src, c.simSel.c_str(), c.simVar.c_str(), c.L.tag.c_str(), c.simScale, (int)c.fact);
    auto it = gTmpl.find(key);
    if (it != gTmpl.end()) return it->second;
    Tmpl& T = gTmpl[key];
    bool ok = c.src == 0 ? LoadTmplF10(gFH, c.simSel.c_str(), c.simVar.c_str(), T, gLog) : LoadTmplG3(gFH, c.simSel.c_str(), c.L, T, gLog);
    if (!ok) P("  !! templates %s failed\n", key.c_str());
    if (c.fact) FactoriseSignal(T, c.L);
    if (c.simScale != 1.) ScaleT(T, c.simScale);
    return T;
}

// F-main derived quantities from (a, kappa, rho, g) and the all-cell window sums
double DerivF(const std::string& n, const Vec& th, const Tmpl& T) {
    const double a = th[0], k = th[1], r = th[2], g = th[3];
    const double W2B = T.wa[0][kWB], W3B = T.wa[1][kWB], W4P = T.wa[2][kWP], W5P = T.wa[3][kWP];
    if (n == "a_eta") return a;
    if (n == "kappa") return k;
    if (n == "rho_d") return r;
    if (n == "g") return g;
    if (n == "Y_bump_F") return kYf * a * (W2B + k * W3B);
    if (n == "Y_peak_F") return kYf * a * (k * W4P + r * W5P);
    if (n == "Y4") return kYf * a * k * W4P;
    if (n == "Y5") return kYf * a * r * W5P;
    if (n == "DR_F") return ((k * W4P + r * W5P) / (W2B + k * W3B)) / ((W4P + W5P) / (W2B + W3B));
    if (n == "BR_A") return kBRdir * r;
    if (n == "sigma_eta_F") return a * kSigEta * kNormEta;
    return NAN;
}
const char* kFq[] = {"a_eta", "kappa", "rho_d", "g", "Y_bump_F", "Y_peak_F", "Y4", "Y5", "DR_F", "BR_A", "sigma_eta_F"};

// smoothing cache per (template key, mode, order, raw, omega)
std::map<std::string, SmSet> gSm;
SmSet& GetSm(const Cfg& c, const Tmpl& T, const DataL& D, int mode, int ord, bool raw, bool omega) {
    const std::string key = Form("%p|%d|%d|%d|%d", (const void*)&T, mode, ord, (int)raw, (int)omega);
    auto it = gSm.find(key);
    if (it != gSm.end()) return it->second;
    SmSet& S = gSm[key];
    SmoothT(T, D.L, mode, ord, raw, S, omega);
    return S;
}

// ---------------------------------------------------------------------------------------------------------------
Res Run(const Cfg& c) {
    TStopwatch sw;
    Res R; R.tag = c.tag; R.fam = c.fam; R.text = c.text;
    DataL& D = GetData(c);
    Tmpl& T = GetTmpl(c);
    const bool omF = c.boF.omega, omA = c.boA.omega;
    SmSet& SF = GetSm(c, T, D, 0, c.ordF, c.rawF, omF);
    SmSet& SA = GetSm(c, T, D, 1, c.ordA, c.rawA, omA);
    // ---- F-main
    Vec thF(4); Mat CF = Zero(4, 4); double kBF = 0;
    Spec SFm; Model MFm; Vec thFfull;
    if (!c.extF) {
        Model M; BkgOpt bo = c.boF; bo.useE1 = (c.useE & 1) != 0;
        Spec S = SpecMainX(c.eo, bo, M);
        if (c.rawF) M.bbBkg = true;
        GBins B;
        MakeGBins(D, T, SF, false, c.useE, c.fLo, c.fHi, B);
        if (c.bkgFixed && c.nomRef) {
            // period test: background of the nominal F-main (bin model) x scale; kappa, g fixed to nominal
            const Vec& t = c.nomRef->thF;   // a kappa rho g
            for (int i = 0; i < B.n; ++i) {
                const int j = B.j[i];
                const Vec& th0 = c.nomRef->thA;  // unused
                (void)th0;
                B.gd[i] = c.nomRef->q.at("_c" + std::to_string(j)).v * (1. + c.nomRef->q.at("_b").v * B.x[i]) * (B.gd[i] + c.nomRef->kBF * B.gc[i]);
                B.gc[i] = 0.; B.gdo[i] = 0.; B.gco[i] = 0.;
            }
            Spec S2; const int ia = S2.Add("a_eta", 0.4), ik = S2.Add("kappa", t[1]), ir = S2.Add("rho_d", 1.3), ig = S2.Add("g", t[3]);
            S2.fixed[ik] = 1; S2.fixed[ig] = 1;
            S2.terms = {{wp5::kS2, ia, -1, -1}, {wp5::kS3, ia, ik, -1}, {wp5::kS4, ia, ik, -1}, {wp5::kS5, ia, ir, -1}, {wp5::kSeta, ia, ig, -1}};
            const int ic = S2.Add("c_bg", 1.0);
            for (int j = 0; j < kNSlice; ++j) S2.iC[j] = ic;
            S2.iB = S2.Add("b", 0.); S2.fixed[S2.iB] = 1;
            S2.iKB = S2.Add("kappa_B", 0.); S2.fixed[S2.iKB] = 1;
            S2.iK[1] = S2.Add("k_E2", 1.0); S2.cons.push_back({S2.iK[1], 1.0, 0.05});
            S = S2; M = Model();
        }
        FitRes F = FitX(B, S, M, D);
        if (!F.ok) { P("  !! [%s] F-main failed\n", c.tag.c_str()); return R; }
        SFm = S; MFm = M; thFfull = F.th;
        for (int k = 0; k < 4; ++k) { thF[k] = F.th[k]; for (int l = 0; l < 4; ++l) CF[k][l] = F.Cstat[k][l]; }
        kBF = M.kBtie ? F.th[1] : (S.iKB >= 0 ? F.th[S.iKB] : 0.);
        R.chiF = F.chi2P; R.ndfF = F.ndf; R.thF = F.th;
        for (int j = 0; j < kNSlice; ++j) if (S.iC[j] >= 0) R.q["_c" + std::to_string(j)].v = F.th[S.iC[j]];
        if (M.iC0 >= 0) for (int j = 0; j < kNSlice; ++j) R.q["_c" + std::to_string(j)].v = F.th[M.iC0] * (1. + F.th[M.iC1] * (D.mpp[j] - 0.38) / 0.05);
        R.q["_b"].v = S.iB >= 0 ? F.th[S.iB] : 0.;
        for (const char* n : kFq) {
            const std::string nm = n;
            auto f = [&](const Vec& th) { Vec t4(th.begin(), th.begin() + 4); return DerivF(nm, t4, T); };
            QV v; v.v = f(F.th);
            const Vec g = wp5::NumGrad(f, F.th, F.C);
            v.stat = std::sqrt(std::max(0., Quad(F.Cstat, g)));
            v.eN = std::sqrt(std::max(0., Quad(F.Cn, g))); v.eCB = std::sqrt(std::max(0., Quad(F.Ccb, g)));
            R.q[nm] = v;
        }
    } else {
        thF = c.thFext; CF = c.CFext; kBF = c.kBext;
        R.thF = thF;
        for (const char* n : kFq) {
            const std::string nm = n;
            auto f = [&](const Vec& th) { return DerivF(nm, th, T); };
            QV v; v.v = f(thF);
            Mat C4 = CF; Vec g(4, 0.);
            for (int k = 0; k < 4; ++k) { const double h = std::max(1e-7, 1e-4 * std::sqrt(std::max(C4[k][k], 0.))); Vec a = thF, b = thF; a[k] += h; b[k] -= h; g[k] = (f(a) - f(b)) / (2 * h); }
            v.stat = std::sqrt(std::max(0., Quad(C4, g)));
            R.q[nm] = v;
        }
    }
    R.CF = CF; R.kBF = kBF;
    if (c.bkgFixed) { R.ok = true; return R; }
    // ---- F-gold
    Spec SGo; Model MGo; Vec thGfull;
    if (c.gold) {
        Model M;
        Spec S = SpecGoldX(c.eo, c.boF, kBF, M);
        if (c.rawF) M.bbBkg = true;
        GBins B;
        const Vec f4 = thF;
        MakeGBins(D, T, SF, true, 2, c.fLo, c.fHi, B, nullptr, [&](int f5, int q, int x) { (void)f5; const size_t i = (size_t)q * T.nb + x; return f4[0] * f4[1] * (T.s5[1][i] + T.s5[2][i]); });
        FitRes F = FitX(B, S, M, D);
        if (F.ok) {
            auto f = [&](const Vec& th) { return kBRdir * th[1] / th[0]; };
            QV v; v.v = f(F.th);
            const Vec g = wp5::NumGrad(f, F.th, F.C);
            v.stat = std::sqrt(std::max(0., Quad(F.Cstat, g)));
            R.q["BR_G"] = v;
            R.chiG = F.chi2P; R.ndfG = F.ndf; R.thG = F.th;
            SGo = S; MGo = M; thGfull = F.th;
        } else P("  !! [%s] F-gold failed\n", c.tag.c_str());
    }
    // ---- stage A (eta fixed to F-main)
    BkgOpt boA = c.boA;
    if (c.aKBfromKappa) { boA.fixKB = true; boA.kBval = thF[1]; }
    Model MA;
    Spec SAp = SpecStageAX(boA, MA);
    if (c.rawA) MA.bbBkg = true;
    auto etaFixT = [](const Tmpl* Tp, const Vec& f) {
        return [Tp, f](int f5, int q, int x) {
            (void)f5; const size_t i = (size_t)q * Tp->nb + x;
            return f[0] * Tp->s[0][i] + f[0] * f[1] * (Tp->s[1][i] + Tp->s[2][i]) + f[0] * f[2] * Tp->s[3][i] + f[0] * f[3] * Tp->s[4][i];
        };
    };
    auto etaFix = [&](const Vec& f) { return etaFixT(&T, f); };
    std::function<bool(int, int, double, double)> keep = nullptr;
    if (c.crk) keep = [](int j, int e, double lo, double hi) { return j == 4 && e == 1 && lo >= 0.400 - 1e-9 && hi <= 0.5225 + 1e-9; };
    auto stageA = [&](const Vec& f, const Vec* start, bool sand, FitRes& F) {
        GBins B;
        MakeGBins(D, T, SA, false, 2, c.aLo, c.aHi, B, keep, etaFix(f));
        F = FitX(B, SAp, MA, D, start, sand);
        return B.n;
    };
    FitRes FA;
    stageA(thF, nullptr, true, FA);
    if (!FA.ok) { P("  !! [%s] stage A failed\n", c.tag.c_str()); return R; }
    R.chiA = FA.chi2P; R.ndfA = FA.ndf; R.thA = FA.th;
    const Vec BA0 = BATable(SA, SAp, MA, FA.th, D);
    // perturbed tables: stage-A parameters and F-main parameters (stage A refitted)
    std::vector<Vec> BAp, BAm; std::vector<double> hA;
    for (int k = 0; k < SAp.np(); ++k) {
        if (SAp.fixed[k] || FA.Cstat[k][k] <= 0) { BAp.push_back(BA0); BAm.push_back(BA0); hA.push_back(0); continue; }
        const double h = 0.05 * std::sqrt(FA.Cstat[k][k]);
        Vec tp = FA.th, tm = FA.th; tp[k] += h; tm[k] -= h;
        BAp.push_back(BATable(SA, SAp, MA, tp, D)); BAm.push_back(BATable(SA, SAp, MA, tm, D)); hA.push_back(h);
    }
    std::vector<Vec> BFp, BFm; Vec hF(4, 0.);
    for (int k = 0; k < 4; ++k) {
        hF[k] = 0.25 * std::sqrt(std::max(CF[k][k], 1e-12));
        Vec tp = thF, tm = thF; tp[k] += hF[k]; tm[k] -= hF[k];
        FitRes Fp, Fm;
        stageA(tp, &FA.th, false, Fp); stageA(tm, &FA.th, false, Fm);
        BFp.push_back(Fp.ok ? BATable(SA, SAp, MA, Fp.th, D) : BA0); BFm.push_back(Fm.ok ? BATable(SA, SAp, MA, Fm.th, D) : BA0);
    }
    // ---- estimators
    struct EN { const char* n; Est f; };
    std::vector<EN> ests = {{"Y_bump", MakeB4All(c.wU)}, {"Y_bump_E2", MakeB4(Cells(1), c.wU)}, {"Y_peak", MakePSB()}, {"a_P", MakeAP()},
                            {"a_B", MakeAB(c.wU)}, {"DR", MakeDR(c.wU)},
                            {"sigma_eta", [](const St& s) { return MakeAP()(s) * kSigEta * kNormEta; }}};
    const int nI = D.nI;
    for (auto& e : ests) {
        St s0{&D.N, &D.CB, &BA0, &T, thF, &D};
        QV v; v.v = e.f(s0);
        // N and CB: exact finite differences over the window items
        std::vector<int> wi; for (int h = 0; h < 10; ++h) for (int w = 0; w < kNW; ++w) wi.push_back(D.IWin(h, w));
        Vec gN(wi.size(), 0.), gC(wi.size(), 0.);
        Vec N = D.N, CB = D.CB;
        for (size_t a = 0; a < wi.size(); ++a) {
            N[wi[a]] += 1.; St s{&N, &D.CB, &BA0, &T, thF, &D}; gN[a] = e.f(s) - v.v; N[wi[a]] -= 1.;
            CB[wi[a]] += 1.; St t{&D.N, &CB, &BA0, &T, thF, &D}; gC[a] = e.f(t) - v.v; CB[wi[a]] -= 1.;
        }
        double vN = 0, vC = 0;
        for (size_t a = 0; a < wi.size(); ++a) {
            if (gN[a] != 0) for (size_t b = 0; b < wi.size(); ++b) if (gN[b] != 0) vN += gN[a] * D.covN[(size_t)wi[a] * nI + wi[b]] * gN[b];
            if (gC[a] != 0) for (size_t b = 0; b < wi.size(); ++b) if (gC[b] != 0) vC += gC[a] * D.covCB[(size_t)wi[a] * nI + wi[b]] * gC[b];
        }
        v.eN = std::sqrt(std::max(0., vN)); v.eCB = std::sqrt(std::max(0., vC));
        Vec gA(SAp.np(), 0.);
        for (int k = 0; k < SAp.np(); ++k) {
            if (hA[k] <= 0) continue;
            St a{&D.N, &D.CB, &BAp[k], &T, thF, &D}, b{&D.N, &D.CB, &BAm[k], &T, thF, &D};
            gA[k] = (e.f(a) - e.f(b)) / (2 * hA[k]);
        }
        v.eA = std::sqrt(std::max(0., Quad(FA.Cstat, gA)));
        Vec gF(4, 0.);
        for (int k = 0; k < 4; ++k) {
            Vec tp = thF, tm = thF; tp[k] += hF[k]; tm[k] -= hF[k];
            St a{&D.N, &D.CB, &BFp[k], &T, tp, &D}, b{&D.N, &D.CB, &BFm[k], &T, tm, &D};
            gF[k] = (e.f(a) - e.f(b)) / (2 * hF[k]);
        }
        v.eF = std::sqrt(std::max(0., Quad(CF, gF)));
        v.stat = std::sqrt(v.eN * v.eN + v.eCB * v.eCB + v.eA * v.eA + v.eF * v.eF);
        R.q[e.n] = v;
    }
    // ---- MC: the same chain with every replica of the simulation (WP5/WP6 procedure: F-main replica r -> F-gold and stage A
    //      with replica r -> estimators with the replica-r window sums; WP4 smoothing redone, WP3 smoothing from the nominal
    //      domain and parameters; fits start at the nominal values, no sandwich)
    const std::vector<Tmpl>& RR = c.simSel == "pre" ? gRepPre : gRep;
    if (c.mc && !RR.empty()) {
        TStopwatch sw2;
        std::map<std::string, Vec> rv;
        auto put = [&](const std::string& n, size_t r, double x) { Vec& v = rv[n]; if (v.empty()) v.assign(RR.size(), NAN); v[r] = x; };
        const bool isNom = c.tag == "nom";
        if (isNom) gNomRepThF.assign(RR.size(), Vec(4, NAN));
        int nok = 0;
        for (size_t r = 0; r < RR.size(); ++r) {
            Tmpl Tr = RR[r];
            if (c.fact) FactoriseSignal(Tr, c.L);
            SmSet SFr, SAr;
            SmoothT(Tr, D.L, 0, c.ordF, c.rawF, SFr, omF);
            SmoothT(Tr, D.L, 1, c.ordA, c.rawA, SAr, omA, &SA);
            Vec thFr = thF;
            if (!c.extF) {
                GBins Br; MakeGBins(D, Tr, SFr, false, c.useE, c.fLo, c.fHi, Br);
                FitRes Fr = FitX(Br, SFm, MFm, D, &thFfull, false);
                if (!Fr.ok) continue;
                for (int k = 0; k < 4; ++k) thFr[k] = Fr.th[k];
                for (const char* n : kFq) put(n, r, DerivF(n, thFr, Tr));
            } else if (c.fixFrep) {
                if (r >= gNomRepThF.size() || !(gNomRepThF[r][0] == gNomRepThF[r][0])) continue;
                thFr = gNomRepThF[r];
            }
            if (isNom) gNomRepThF[r] = thFr;
            if (c.gold && !thGfull.empty()) {
                GBins Bg;
                MakeGBins(D, Tr, SFr, true, 2, c.fLo, c.fHi, Bg, nullptr, [&Tr, thFr](int f5, int q, int x) { (void)f5; const size_t i = (size_t)q * Tr.nb + x; return thFr[0] * thFr[1] * (Tr.s5[1][i] + Tr.s5[2][i]); });
                FitRes Fg = FitX(Bg, SGo, MGo, D, &thGfull, false);
                if (Fg.ok) put("BR_G", r, kBRdir * Fg.th[1] / Fg.th[0]);
            }
            GBins Ba;
            MakeGBins(D, Tr, SAr, false, 2, c.aLo, c.aHi, Ba, keep, etaFixT(&Tr, thFr));
            FitRes Far = FitX(Ba, SAp, MA, D, &FA.th, false);
            if (!Far.ok) continue;
            const Vec BAr = BATable(SAr, SAp, MA, Far.th, D);
            St s{&D.N, &D.CB, &BAr, &Tr, thFr, &D};
            for (auto& e : ests) put(e.n, r, e.f(s));
            ++nok;
        }
        for (auto& kv : rv) {
            R.q[kv.first].rep = kv.second;
            Vec x; for (double y : kv.second) if (y == y) x.push_back(y);
            const int n = (int)x.size();
            if (n < 3) continue;
            double m = 0; for (double y : x) m += y; m /= n;
            double v2 = 0; for (double y : x) v2 += (y - m) * (y - m);
            R.q[kv.first].mc = std::sqrt(v2 / (n - 1));
            R.q[kv.first].nrep = n;
        }
        P("  [%s] replica chain: %d of %zu replicas (%s__nom), %.1f s\n", c.tag.c_str(), nok, RR.size(), c.simSel.c_str(), sw2.RealTime());
    }
    R.ok = true;
    P("  [%s] %.1f s\n", c.tag.c_str(), sw.RealTime());
    return R;
}

void PrintRes(const Res& r) {
    if (!r.ok) { P("%-16s FAILED\n", r.tag.c_str()); return; }
    auto g = [&](const char* n) { auto it = r.q.find(n); return it == r.q.end() ? QV() : it->second; };
    P("%-16s %-4s Yb %8.1f+-%6.1f YbE2 %7.1f+-%6.1f Yp %7.1f+-%5.1f aP %.4f aB %.4f DR %.3f | a %.4f+-%.4f k %.3f r %.3f g %.3f Y4 %.0f Y5 %.0f BRA %.3e BRG %.3e DRF %.3f | chi2 F %.1f/%d A %.1f/%d G %.1f/%d\n",
      r.tag.c_str(), r.fam.c_str(), g("Y_bump").v, g("Y_bump").stat, g("Y_bump_E2").v, g("Y_bump_E2").stat, g("Y_peak").v, g("Y_peak").stat, g("a_P").v, g("a_B").v,
      g("DR").v, g("a_eta").v, g("a_eta").stat, g("kappa").v, g("rho_d").v, g("g").v, g("Y4").v, g("Y5").v, g("BR_A").v, g("BR_G").v, g("DR_F").v,
      r.chiF, r.ndfF, r.chiA, r.ndfA, r.chiG, r.ndfG);
}
void WriteRes(const Res& r) {
    if (!gRaw) return;
    if (!r.ok) { fprintf(gRaw, "cfg %s %s FAILED\n", r.tag.c_str(), r.fam.c_str()); return; }
    fprintf(gRaw, "cfg %s %s chiF %.3f %d chiA %.3f %d chiG %.3f %d | %s\n", r.tag.c_str(), r.fam.c_str(), r.chiF, r.ndfF, r.chiA, r.ndfA, r.chiG, r.ndfG, r.text.c_str());
    for (const auto& kv : r.q) {
        if (kv.first[0] == '_') continue;
        const QV& v = kv.second;
        fprintf(gRaw, "q %s %s %.10g %.6g N %.6g CB %.6g A %.6g F %.6g mc %.6g %d\n", r.tag.c_str(), kv.first.c_str(), v.v, v.stat, v.eN, v.eCB, v.eA, v.eF, v.mc, v.nrep);
        if (gRepF && !v.rep.empty()) { fprintf(gRepF, "rep %s %s", r.tag.c_str(), kv.first.c_str()); for (double x : v.rep) fprintf(gRepF, " %.8g", x); fprintf(gRepF, "\n"); }
    }
    fflush(gRaw);
}
}  // namespace wp7

void wp7_syst(int only = 0) {
    using namespace alldefs;
    using namespace wp7;
    TStopwatch total;
    gSystem->mkdir(kOut, kTRUE);
    const char* sfx = only == 1 ? "_check" : "";   // the nominal-only check never overwrites the full outputs
    gLog = fopen(Form("%s/wp7_syst%s.log", kOut, sfx), "w");
    gRaw = fopen(Form("%s/syst_raw%s.txt", kOut, sfx), "w");
    gRepF = fopen(Form("%s/syst_rep%s.txt", kOut, sfx), "w");
    fprintf(gRepF, "# rep <tag> <quantity> value of replica 0..99 (nan: fit failed); same replica index = same simulation replica\n");
    fprintf(gRaw, "# wp7_syst.C (SPEC report_all, WP7), 2026-11-14. S0 nominal, isBest >= 0, L = %.1f nb^-1; yields raw rows x 64/63.\n", kL);
    fprintf(gRaw, "# cfg <tag> <family> chiF chiA chiG | text\n# q <tag> <quantity> value stat N CB A(stage A) F(F-main) mc(replica RMS, nan if not run) nrep\n");
    if (!LoadRows(gRows)) { P("rows missing\n"); return; }
    P("rows loaded: %zu\n", gRows.n);
    gFH = TFile::Open(kHistFile);
    {
        Cfg c0; Tmpl& T0 = GetTmpl(c0);
        if (only != 1 || true) { if (!LoadTmplF10Rep(gFH, "S0", T0, 100, gRep)) { P("replicas missing\n"); gRep.clear(); } }
        Cfg cp; cp.simSel = "pre"; Tmpl& Tp = GetTmpl(cp);
        if (!LoadTmplF10Rep(gFH, "pre", Tp, 100, gRepPre)) { P("pre replicas missing\n"); gRepPre.clear(); }
        P("replicas loaded: %zu (S0), %zu (pre)\n", gRep.size(), gRepPre.size());
    }

    // ================================================================ nominal (and checks against WP5 / WP6)
    Cfg nom; nom.tag = "nom"; nom.fam = "nom"; nom.mc = true; nom.text = "nominal: S0, CB spread rule, F10 templates S0__nom, F-main E2 0.40-0.70, stage A CR-U 0.575-0.700";
    Res rNom = Run(nom);
    PrintRes(rNom); WriteRes(rNom);
    {
        P("CHECK vs WP5 results.txt / WP6 results.txt:\n");
        P("  a_eta %.8f (WP5 0.38228397), kappa %.8f (0.83233498), rho_d %.7f (1.3152029), g %.7f (1.7167655)\n", rNom.q["a_eta"].v, rNom.q["kappa"].v, rNom.q["rho_d"].v, rNom.q["g"].v);
        P("  Y_bump_F %.3f (19426.131), Y_peak_F %.3f (5044.7742), Y4 %.3f (1847.0205), Y5 %.3f (3197.7536), DR_F %.6f (1.1689618), BR_A %.6e (3.5247438e-4)\n",
          rNom.q["Y_bump_F"].v, rNom.q["Y_peak_F"].v, rNom.q["Y4"].v, rNom.q["Y5"].v, rNom.q["DR_F"].v, rNom.q["BR_A"].v);
        P("  stat: a_eta %.6f (0.016925), Y5 %.2f (290.47), BR_A %.4e (3.083e-5)\n", rNom.q["a_eta"].stat, rNom.q["Y5"].stat, rNom.q["BR_A"].stat);
        P("  BR_G %.6e (WP5 2.6546e-4, stat 0.368e-4: %.4e)\n", rNom.q["BR_G"].v, rNom.q["BR_G"].stat);
        P("  B4 E2 %.3f +- %.3f (WP6 7822.17 +- 456.215), B4 all %.3f +- %.3f (16389.3 +- 1129.04), PSB %.3f +- %.3f (5087.97 +- 431.626)\n",
          rNom.q["Y_bump_E2"].v, rNom.q["Y_bump_E2"].stat, rNom.q["Y_bump"].v, rNom.q["Y_bump"].stat, rNom.q["Y_peak"].v, rNom.q["Y_peak"].stat);
        P("  a_P %.6f +- %.6f (0.418249 +- 0.0354812), a_B %.6f (0.311677), DR %.6f +- %.6f (1.39743 +- 0.152024)\n", rNom.q["a_P"].v, rNom.q["a_P"].stat,
          rNom.q["a_B"].v, rNom.q["DR"].v, rNom.q["DR"].stat);
        P("  parts B4 E2: N %.2f CB %.2f A %.2f F %.2f (WP6 250.378 87.1879 368.265 47.136)\n", rNom.q["Y_bump_E2"].eN, rNom.q["Y_bump_E2"].eCB, rNom.q["Y_bump_E2"].eA, rNom.q["Y_bump_E2"].eF);
        P("  MC (replica RMS): a_eta %.6f (WP5 0.039629), Y_bump_F %.1f (2274.1), Y5 %.2f (570.66), BR_G %.4e (rms 6.6354e-5), B4 E2 %.1f (WP6 1919.6), B4 all %.1f (4109.78), PSB %.2f (245.122), a_P %.5f (0.0248145), DR %.4f (0.533303)\n",
          rNom.q["a_eta"].mc, rNom.q["Y_bump_F"].mc, rNom.q["Y5"].mc, rNom.q["BR_G"].mc, rNom.q["Y_bump_E2"].mc, rNom.q["Y_bump"].mc, rNom.q["Y_peak"].mc, rNom.q["a_P"].mc, rNom.q["DR"].mc);
        P("  stage A: ");
        for (double x : rNom.thA) P(" %.5f", x);
        P("  (WP6: 0.9816 0.9844 0.9505 1.0473 0.9145 -0.0639 2.932 0.955)\n");
    }
    // F-main of the nominal for the configurations that keep it fixed
    if (only == 1) { fclose(gRaw); gRaw = nullptr; fclose(gRepF); gRepF = nullptr; fclose(gLog); gLog = nullptr; return; }
    auto fixF = [&](Cfg& c) { c.extF = true; c.thFext = Vec(rNom.thF.begin(), rNom.thF.begin() + 4); c.kBext = rNom.kBF; c.CFext = rNom.CF; c.gold = false; c.fixFrep = true; };
    std::vector<Res> all; all.push_back(rNom);
    auto run = [&](Cfg c) {
        // the replica chain only where the variant uses the S0__nom simulation sample in the standard layout
        c.mc = c.src == 0 && (c.simSel == "S0" || c.simSel == "pre") && c.simVar == "nom" && c.simScale == 1. && c.L.tag == "std" && !c.bkgFixed && c.tag != "S4_GLShyb" &&
               c.fam != "S1info" && c.fam != "S2info" && c.fam != "test" && c.fam != "ref";
        Res r = Run(c); PrintRes(r); WriteRes(r); all.push_back(r); return r; };

    // reference for the G3-built layouts: the standard layout from G3
    Cfg g3 = nom; g3.tag = "nomG3"; g3.fam = "ref"; g3.src = 1; g3.text = "nominal with the simulation summed from the G3 layout (reference of the G3-based variants)";
    run(g3);

    // ================================================================ S1 CB
    { Cfg c = nom; c.tag = "S1_kE2fix1"; c.fam = "S1"; c.boF.kE2fix = 1.0; c.boA.kE2fix = 1.0; c.text = "k_E2 fixed to 1 (F-main, F-gold, stage A)"; run(c); }
    { Cfg c = nom; c.tag = "S1_k1032"; c.fam = "S1"; c.C.k = 1.032; c.boF.kE2fix = 1.0; c.boA.kE2fix = 1.0; c.text = "CB x 1.032 everywhere (k from data mixing), k_E2 fixed"; run(c); }
    { Cfg c = nom; c.tag = "S1_LSsum"; c.fam = "S1"; c.C.rule = 1; c.text = "CB = N++ + N-- (like-sign arithmetic sum) instead of 2 sqrt(N++ N--)"; run(c); }
    { Cfg c = nom; c.tag = "S1_CB20"; c.fam = "S1"; c.C.reb = kCBreb20; c.text = "CB in 20 MeV cells"; run(c); }
    { Cfg c = nom; c.tag = "S1_geo"; c.fam = "S1"; c.C.rule = 2; c.text = "CB per-cell geometric mean (10 MeV cells of every histogram)"; run(c); }
    { Cfg c = nom; c.tag = "S1_E1kept"; c.fam = "S1info"; c.useE = 3; c.eo.kappaFree = true; c.text = "INFO: F-main as specified, E1 + E2 with k_E1 free (rejected: k_E1 outside 0.8-3)"; run(c); }

    // ================================================================ S2 CR definition (stage A and the side-band estimators; F-main fixed)
    { Cfg c = nom; fixF(c); c.tag = "S2_U560"; c.fam = "S2"; c.aLo = 0.560; c.wU = kWU560; c.text = "CR-U 0.560-0.700 (stage A and B4 U)"; run(c); }
    { Cfg c = nom; fixF(c); c.tag = "S2_U565"; c.fam = "S2"; c.aLo = 0.565; c.wU = kWU565; c.text = "CR-U 0.565-0.700"; run(c); }
    { Cfg c = nom; fixF(c); c.tag = "S2_Uhi650"; c.fam = "S2"; c.aHi = 0.650; c.wU = kWU650; c.text = "CR-U 0.575-0.650"; run(c); }
    Layout L75 = StdLayout(); L75.tag = "r075"; L75.e = EdgesF10like(0.750);
    Layout L85 = StdLayout(); L85.tag = "r085"; L85.e = EdgesF10like(0.850);
    { Cfg c = g3; fixF(c); c.tag = "S2_Uhi750"; c.fam = "S2"; c.L = L75; c.aHi = 0.750; c.wU = kWU750; c.text = "CR-U 0.575-0.750 (G3 layout to 0.75)"; run(c); }
    { Cfg c = nom; fixF(c); c.tag = "S2_CRK043"; c.fam = "S2info"; c.crk = true; c.text = "INFO: stage A with CR-K (S5 E2, 0.400-0.5225) added (CR-K rejected by A1/A2)"; run(c); }
    { Cfg c = g3; fixF(c); c.tag = "S2_CRK044"; c.fam = "S2info"; c.L.tag = "s5_044"; c.L.se[3] = 0.44; c.crk = true; c.text = "INFO: S4/S5 edge at 0.44 and CR-K added (G3)"; run(c); }
    { Cfg c = g3; fixF(c); c.tag = "S2_CRK0425"; c.fam = "S2info"; c.L.tag = "s5_0425"; c.L.se[3] = 0.425; c.crk = true; c.text = "INFO: S4/S5 edge at 0.425 and CR-K added (G3)"; run(c); }
    { Cfg c = g3; fixF(c); c.tag = "S2_K0veto"; c.fam = "S2"; c.L.tag = "veto"; c.L.veto = true; c.text = "K0S veto M(pipi) 0.485-0.507 (data and simulation, G3)"; run(c); }

    // ================================================================ S3 background model (F-main and stage A)
    { Cfg c = nom; c.tag = "S3_quad"; c.fam = "S3"; c.boF.quad = c.boA.quad = true; c.text = "+ quadratic term q x^2 (F-main and stage A)"; run(c); }
    { Cfg c = nom; c.tag = "S3_tiltj"; c.fam = "S3"; c.boF.tiltj = c.boA.tiltj = true; c.text = "per-slice tilts b_j"; run(c); }
    { Cfg c = nom; c.tag = "S3_expb"; c.fam = "S3"; c.boF.expb = c.boA.expb = true; c.text = "tilt exp(b x) instead of 1 + b x"; run(c); }
    { Cfg c = nom; c.tag = "S3_commonC"; c.fam = "S3"; c.boF.commonC = c.boA.commonC = true; c.text = "common scale c0 (1 + c1 y), y = (<M(pipi)>_slice - 0.38)/0.05"; run(c); }
    { Cfg c = nom; c.tag = "S3_omega"; c.fam = "S3"; c.boF.omega = c.boA.omega = true; c.text = "omega/rho0-pion background with its own free scale c_om"; run(c); }
    { Cfg c = g3; c.tag = "S3_omega085"; c.fam = "S3"; c.L = L85; c.boF.omega = c.boA.omega = true; c.fHi = 0.85; c.aHi = 0.85; c.wU = kWU850; c.text = "omega split and range to 0.85 (F-main 0.40-0.85, CR-U 0.575-0.85; G3)"; run(c); }
    { Cfg c = nom; c.tag = "S3_kBtie"; c.fam = "S3"; c.boF.kBtie = true; c.aKBfromKappa = true; c.text = "kappa_B = kappa (F-main); stage A kappa_B fixed to the F-main kappa"; run(c); }
    { Cfg c = nom; c.tag = "S3_D11data"; c.fam = "S3"; c.boF.tiltj = c.boA.tiltj = true; c.boF.quad = c.boA.quad = true; c.text = "WP3 data-only D11 choice: per-slice tilts + quadratic"; run(c); }
    { Cfg c = nom; c.tag = "S3_D11dataOm"; c.fam = "S3"; c.boF.tiltj = c.boA.tiltj = true; c.boF.quad = c.boA.quad = true; c.boF.omega = c.boA.omega = true; c.text = "per-slice tilts + quadratic + omega split"; run(c); }

    // ================================================================ S4 background template
    { Cfg c = nom; c.tag = "S4_pol4"; c.fam = "S4"; c.ordF = c.ordA = 4; c.text = "exp(pol4) smoothing of G_D, G_C (F-main and stage A)"; run(c); }
    { Cfg c = nom; c.tag = "S4_raw"; c.fam = "S4"; c.rawF = c.rawA = true; c.text = "raw G_D, G_C with BB-lite (their variance in the weights)"; run(c); }
    for (const char* v : {"ghostA", "modeB", "wdb"}) { Cfg c = nom; c.tag = std::string("S4_") + v; c.fam = "S4"; c.simVar = v; c.text = std::string("simulation variant ") + v; run(c); }
    {
        // template-noise-aware estimator (v5 hybrid GLS, minimised; parameters from v5_quad2.log): F-main replaced
        Cfg c = nom; c.tag = "S4_GLShyb"; c.fam = "S4"; c.extF = true; c.thFext = {0.5126, 0.856, 1.118, 0.465}; c.kBext = 2.871; c.CFext = rNom.CF; c.gold = true;
        c.text = "F-main with the hybrid GLS (data + BB-lite signal + replica background covariance, v5): a 0.5126, kappa 0.856, rho 1.118, g 0.465, kappa_B 2.871";
        run(c);
    }

    // ================================================================ S5 range, binning
    { Cfg c = nom; c.tag = "S5_r4265"; c.fam = "S5"; c.fLo = 0.42; c.fHi = 0.65; c.text = "F-main and F-gold on 0.42-0.65"; run(c); }
    { Cfg c = nom; c.tag = "S5_r4570"; c.fam = "S5"; c.fLo = 0.45; c.text = "F-main and F-gold on 0.45-0.70"; run(c); }
    { Cfg c = g3; c.tag = "S5_r4075"; c.fam = "S5"; c.L = L75; c.fHi = 0.75; c.text = "F-main and F-gold on 0.40-0.75 (G3)"; run(c); }
    { Cfg c = g3; c.tag = "S5_bin25"; c.fam = "S5"; c.L.tag = "b25"; c.L.e = EdgesF10like(0.700, 0.0025, 0.0025, 0.0025); c.text = "2.5 MeV bins everywhere (G3)"; run(c); }
    { Cfg c = g3; c.tag = "S5_sliceUp"; c.fam = "S5"; c.L.tag = "slup"; for (int q = 0; q < 4; ++q) c.L.se[q] += 0.005; c.text = "slice edges + 5 MeV (G3)"; run(c); }
    { Cfg c = g3; c.tag = "S5_sliceDn"; c.fam = "S5"; c.L.tag = "sldn"; for (int q = 0; q < 4; ++q) c.L.se[q] -= 0.005; c.text = "slice edges - 5 MeV (G3)"; run(c); }
    { Cfg c = g3; c.tag = "S5_mee08"; c.fam = "S5"; c.L.tag = "mee08"; c.L.eSplit = 0.008; c.text = "M(ee) split at 8 MeV (G3)"; run(c); }
    { Cfg c = g3; c.tag = "S5_mee12"; c.fam = "S5"; c.L.tag = "mee12"; c.L.eSplit = 0.012; c.text = "M(ee) split at 12 MeV (G3)"; run(c); }

    // ================================================================ S6 mass transfer
    for (const char* v : {"sUp", "sDn", "sK0", "r089", "r110"}) { Cfg c = nom; c.tag = std::string("S6_") + v; c.fam = "S6"; c.simVar = v; c.text = std::string("simulation variant ") + v; run(c); }
    { Cfg c = g3; c.tag = "S6_pipiUp"; c.fam = "S6"; c.L.tag = "pipiup"; c.L.simEdgeRel = -0.0007 / 0.9952; c.text = "simulated M(pipi) x (0.9952 + 0.0007)/0.9952 (slice edges shared linearly in the 5 MeV G3 bins)"; run(c); }
    { Cfg c = g3; c.tag = "S6_pipiDn"; c.fam = "S6"; c.L.tag = "pipidn"; c.L.simEdgeRel = +0.0007 / 0.9952; c.text = "simulated M(pipi) x (0.9952 - 0.0007)/0.9952"; run(c); }

    // ================================================================ S7 physics weights
    for (const char* v : {"V2", "bump1aY", "bumpPS", "convOff", "GS"}) { Cfg c = nom; c.tag = std::string("S7_") + v; c.fam = "S7"; c.simVar = v; c.text = std::string("physics weights ") + v; run(c); }
    { Cfg c = nom; c.tag = "S7_wtsOff"; c.fam = "S7info"; c.simVar = "wtsOff"; c.text = "INFO: physics weights off (shown only)"; run(c); }

    // ================================================================ S8 G_eta
    { Cfg c = nom; c.tag = "S8_g015"; c.fam = "S8"; c.eo.gSig = 0.15; c.text = "g = 1 +- 0.15"; run(c); }
    { Cfg c = nom; c.tag = "S8_g050"; c.fam = "S8"; c.eo.gSig = 0.50; c.text = "g = 1 +- 0.50"; run(c); }
    { Cfg c = nom; c.tag = "S8_gfree"; c.fam = "S8"; c.eo.gFree = true; c.text = "g free"; run(c); }

    // ================================================================ S9 signal templates
    { Cfg c = nom; c.tag = "S9_fact"; c.fam = "S9"; c.fact = true; c.text = "factorised signal templates f(ee) g(slice) h(M4 | slice)"; run(c); }

    // ================================================================ S11 selection
    { Cfg c = nom; c.tag = "S11_pre"; c.fam = "S11"; c.S.sel = kSelPre; c.simSel = "pre"; c.text = "selection pre (no pion-vertex requirement), pre templates"; run(c); }
    { Cfg c = nom; c.tag = "S11_dz40"; c.fam = "S11"; c.S.sel = kSelDz40; c.simSel = "dz40"; c.text = "pre + dz_pi < 40 mm instead of the pion-vertex flags"; run(c); }
    { Cfg c = nom; c.tag = "S11_vz150"; c.fam = "S11"; c.S.sel = kSelVz150; c.simSel = "vz150"; c.text = "S0 with vz > -150 mm"; run(c); }
    { Cfg c = nom; c.tag = "S11_vz170"; c.fam = "S11"; c.S.sel = kSelVz170; c.simSel = "vz170"; c.text = "S0 with vz > -170 mm"; run(c); }
    { Cfg c = nom; c.tag = "S11_PEAK"; c.fam = "S11info"; c.S.sel = kSelPeak; c.simSel = "PEAK"; c.text = "INFO: PEAK = pre + dz_pi < 40 + MM4 < 2.85 (D4 variant, not counted)"; run(c); }

    // ================================================================ consistency tests: folds, isBest == 1, cor_ok fallback
    double nAll = 0;
    { DataL& D = GetData(nom); nAll = D.Ntot; }
    for (int f = 0; f < 2; ++f) {
        Cfg c = nom; c.S.sub = f == 0 ? kSubFold0 : kSubFold1; c.tag = f == 0 ? "T_fold0" : "T_fold1"; c.fam = "test";
        DataL& D = GetData(c);
        c.simScale = D.Ntot / nAll;
        c.text = Form("data fold %d, templates x %.5f (row fraction)", f, c.simScale);
        run(c);
    }
    { Cfg c = nom; c.tag = "T_best"; c.fam = "test"; c.S.best = true; c.text = "data isBest == 1 (templates isBest >= 0, D2)"; run(c); }
    { Cfg c = nom; c.tag = "T_corfb"; c.fam = "test"; c.S.corfb = true; c.text = "data with the cor_ok fallback rows (templates nominal)"; run(c); }

    // ================================================================ periods (WP2 stability test): background fixed
    {
        const double Lp[4] = {2809.453, 372.610, 1078.022, 1819.110};   // analysed L per period (WP2 Table N1), sum 6079.195
        double Ls = 0; for (double x : Lp) Ls += x;
        FILE* pf = fopen(Form("%s/periods_raw.txt", kOut), "w");
        fprintf(pf, "# period L f  | F-main type (bkg shape and kappa, g fixed to nominal; a_eta, rho_d, scale free): a_eta stat DR_F stat | side-band type (B_A shape and F-main fixed to nominal): a_P stat(N+CB) a_B stat Y_peak stat Y_bump stat DR stat\n");
        for (int p = 0; p < 4; ++p) {
            Cfg c = nom; c.S.sub = kSubPerI + p; c.simScale = Lp[p] / Ls;
            // (A) F-main type
            Cfg ca = c; ca.tag = Form("T_per%s_F", PeriodName(p)); ca.fam = "period"; ca.bkgFixed = true; ca.nomRef = &rNom; ca.gold = false;
            ca.text = Form("period %s, F-main with the nominal background shape x one scale; kappa, g fixed", PeriodName(p));
            Res ra = Run(ca); PrintRes(ra); WriteRes(ra);
            // (B) side-band type: nominal stage A (B_A shape), nominal F-main; only the period data change
            DataL& D = GetData(c);
            Tmpl& T = GetTmpl(c);
            DataL& D0 = GetData(nom);
            Tmpl& T0 = GetTmpl(nom);
            SmSet& SA0 = GetSm(nom, T0, D0, 1, 3, false, false);
            Model MA; Spec SAp = SpecStageAX(nom.boA, MA);
            Vec BA0 = BATable(SA0, SAp, MA, rNom.thA, D0);
            for (auto& x : BA0) x *= c.simScale;   // the nominal background per unit L (only the shape matters for R; R' mixes it with the tails)
            const Vec thF(rNom.thF.begin(), rNom.thF.begin() + 4);
            fprintf(pf, "%s %.3f %.6f %.6g %.6g %.6g %.6g", PeriodName(p), Lp[p], c.simScale, ra.q["a_eta"].v, ra.q["a_eta"].stat, ra.q["DR_F"].v, ra.q["DR_F"].stat);
            struct EN { const char* n; Est f; };
            std::vector<EN> ests = {{"a_P", MakeAP()}, {"a_B", MakeAB(kWU575)}, {"Y_peak", MakePSB()}, {"Y_bump", MakeB4All(kWU575)}, {"DR", MakeDR(kWU575)}};
            for (auto& e : ests) {
                St s0{&D.N, &D.CB, &BA0, &T, thF, &D};
                const double v = e.f(s0);
                std::vector<int> wi; for (int h = 0; h < 10; ++h) for (int w = 0; w < kNW; ++w) wi.push_back(D.IWin(h, w));
                Vec gN(wi.size()), gC(wi.size());
                Vec N = D.N, CB = D.CB;
                for (size_t a = 0; a < wi.size(); ++a) {
                    N[wi[a]] += 1; St s{&N, &D.CB, &BA0, &T, thF, &D}; gN[a] = e.f(s) - v; N[wi[a]] -= 1;
                    CB[wi[a]] += 1; St t{&D.N, &CB, &BA0, &T, thF, &D}; gC[a] = e.f(t) - v; CB[wi[a]] -= 1;
                }
                double vv = 0;
                for (size_t a = 0; a < wi.size(); ++a) for (size_t b = 0; b < wi.size(); ++b)
                    vv += gN[a] * D.covN[(size_t)wi[a] * D.nI + wi[b]] * gN[b] + gC[a] * D.covCB[(size_t)wi[a] * D.nI + wi[b]] * gC[b];
                fprintf(pf, " %.6g %.6g", v, std::sqrt(std::max(0., vv)));
                P("  period %s %s = %.4g +- %.3g (N + CB)\n", PeriodName(p), e.n, v, std::sqrt(std::max(0., vv)));
            }
            fprintf(pf, "\n");
        }
        fclose(pf);
    }

    P("\nwp7_syst: %zu configurations, total %.0f s\n", all.size(), total.RealTime());
    fclose(gRaw); gRaw = nullptr;
    fclose(gRepF); gRepF = nullptr;
    fclose(gLog); gLog = nullptr;
}
