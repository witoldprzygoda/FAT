// ============================================================================
// sim_tofmom.C  --  effect of the expected flight time (tof_mom) on beta_new
// ============================================================================
// The simulation PAT (branch pp45_sim_pippimepem, commit 3b8f529, user/htimecut.cc)
// computes the expected flight time in two segments,
//     tof_mom = track_length / C * ( f / beta(p_corr) + (1 - f) / beta(p) ),
//     f = 900 / tracklength   (only if tracklength > 900 mm and p, p_corr > 0),
// the data PAT (tag prod-feb22-betanew-2026-10-28) in one segment,
//     tof_mom = track_length / ( beta(p) C ),   track_length = beta_DST C tof.
// HTimeCut::calcTof (both branches identical):
//     reference i = e- (sort index 0, the only one in PipPimEpEm),
//     tof_new_j = (tm_i + tm_j)/2 - (t_i - t_j)/2                (j != i)
//     tof_new_i = sum_{j != i} [ (tm_i + tm_j)/2 + (t_i - t_j)/2 ] / (size - 1)
//     beta_new  = 1e6 tracklength / (C tof_new)
//     chi2_red  = sum_j (tof_new_ij - tm_j)^2 / resolution_j^2 / (size - 1)
// PID (HGraphCut): TCutG <x>_cut_3sig (pions), <x>_cut_5sig (leptons) in (p, beta_new);
// data cut file pp45_pid_cuts_exp_betanew_v1.root, simulation pp45_pid_cuts_sim_betanew_v6.root,
// both extracted with git show from the PAT repository into
// /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/tofmom/.
//
// Steps
//  1. reproduce the stored tof_mom (data: one segment, sim: two segments), tof_new,
//     beta_new and hyp_best_chi2_red row by row;
//  2. data: two-segment tof_mom -> change of tof_new, beta_new, chi2 for the four
//     particles; delta beta_new vs p; measured beta_new resolution per p bin;
//     tracks crossing the PID polygons (data cut file) and n-sigma bands
//     (3 sigma pions, 5 sigma leptons, measured mu and sigma); mean-shift-removed
//     variant (polygons refitted after a reprocessing); effect on M4 and on isBest;
//  3. simulation: the same quantities for the reverse change (two -> one segment).
//
// Convention: delta = (two segments) - (one segment) on both sides.
//
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_tofmom", kTRUE)' 'plots/eta_sim/sim_tofmom.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TCutG.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

namespace tofmom {

const char* kDataFile = "/hdd3/przygoda/HADES/files/exp/gen4/lep_034.root";
const char* kSimFile  = "/hdd3/przygoda/HADES/files/sim/ppee_std_all.root";
const char* kCutExp   = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/tofmom/pp45_pid_cuts_exp_betanew_v1.root";
const char* kCutSim   = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/tofmom/pp45_pid_cuts_sim_betanew_v6.root";
const char* kFig      = "plots/eta_sim/figs/";
const char* kTab      = "plots/eta_sim/tables/";

const double kC   = 299792458.;   // as HTimeCut (m/s; tof in ns, tracklength in mm)
const double kMpi = 139.56995, kMe = 0.51099906;   // MeV/c2, as HTimeCut
const double kL1  = 900.;         // kElossSegmentLength [mm]
const int    kNP  = 4;
const char*  kPart[kNP]  = {"pip", "pim", "ep", "em"};
const char*  kLabel[kNP] = {"#pi^{+}", "#pi^{-}", "e^{+}", "e^{-}"};
const char*  kCorr[kNP]  = {"p_corr_pip", "p_corr_pim", "p_corr_ep", "p_corr_em"};
const double kMass[kNP]  = {kMpi, kMpi, kMe, kMe};
const double kNsig[kNP]  = {3., 3., 5., 5.};
const int    kRef = 3;            // e- is the time reference (sort index 0)

// momentum bins of the resolution measurement [MeV/c]
const std::vector<double> kPbPi = {60, 90, 120, 150, 180, 220, 260, 300, 350, 400, 500, 600, 800, 1000, 1300, 1700, 2500};
const std::vector<double> kPbE  = {0, 60, 100, 150, 200, 250, 300, 400, 500, 700, 1000, 1900};
const std::vector<double>& Pb(int k) { return k < 2 ? kPbPi : kPbE; }
int PBin(int k, double p) {
    const auto& b = Pb(k);
    if (p < b.front() || p >= b.back()) return -1;
    return int(std::upper_bound(b.begin(), b.end(), p) - b.begin()) - 1;
}

struct PIn { float p, pc, beta, tof, tm, tn, bn, tl, L, res, th, ph; };
struct RIn { PIn x[kNP]; float isBest, run, ev, trig, vz, chi2; };

struct TC { double tn[kNP], bn[kNP], chi2; };

double BetaP(double p, double m) { return p / std::sqrt(p * p + m * m); }
double Tm1(const PIn& x, double m) { return x.tl / (kC * BetaP(x.p, m)); }
double Tm2(const PIn& x, double m) {
    double t = Tm1(x, m);
    if (x.p > 0. && x.pc > 0. && x.L > kL1) {
        const double f = kL1 / x.L;
        t = x.tl / kC * (f / BetaP(x.pc, m) + (1. - f) / BetaP(x.p, m));
    }
    return t;
}
// HTimeCut::calcTof for the PipPimEpEm hypothesis (reference = e-)
TC TimeCut(const double t[kNP], const double tm[kNP], const double L[kNP], const double res[kNP]) {
    TC o; double arr[kNP] = {0, 0, 0, 0};
    const int i = kRef;
    for (int j = 0; j < kNP; ++j) {
        if (j == i) continue;
        const double dt = (t[i] - t[j]) / 2., mean = (tm[i] + tm[j]) / 2.;
        arr[j] = mean - dt;
        arr[i] += (mean + dt) / double(kNP - 1);
    }
    double chi2 = 0;
    for (int j = 0; j < kNP; ++j) chi2 += std::pow(arr[j] - tm[j], 2) / std::pow(res[j], 2);
    o.chi2 = chi2 / double(kNP - 1);
    for (int j = 0; j < kNP; ++j) { o.tn[j] = arr[j]; o.bn[j] = 1e6 * L[j] / (kC * arr[j]); }
    return o;
}

// one row after processing
struct Row {
    bool sel, best;
    long long run, ev;
    double vz;
    double p[kNP], bn0[kNP];      // stored p, stored beta_new
    double dtm[kNP];              // tm(two) - tm(one), formula
    double db[kNP];               // beta_new(two) - beta_new(one)
    double chi2_0, chi2_1, chi2_2;   // stored, one-segment, two-segment (recomputed)
    double m4;                    // GeV/c2, p_corr
    bool in0[kNP], inS[kNP], inM[kNP], inR[kNP], inRM[kNP];   // polygon: stored, shifted, mirror, shifted - mean, mirror of that
    bool newBest = false;         // lowest two-segment chi2 of its event
    double isb;
};

struct Verify {
    long long n = 0;
    double maxd_tm1[kNP] = {0}, maxd_tm2[kNP] = {0};   // |stored - formula|
    double maxd_tn = 0, maxd_bn = 0, maxd_chi2rel = 0;
    long long nbad_tm = 0;    // |diff| > 1e-4 ns with the formula of that production
    long long nbad_bn = 0;    // |diff| > 1e-5
    long long nbad_chi2 = 0;  // relative diff > 1e-3
    long long nin0 = 0, nin_rows = 0; // particles / rows inside the polygons of that production
    long long ntr = 0;
    TH1D* hlog[2];            // log10 |stored - formula| (correct formula, other formula), all particles
    TH1D* hlogbn;             // log10 |stored beta_new - recomputed|
};

TCutG* GetCut(TFile* f, const char* n) {
    TCutG* c = (TCutG*)f->Get(n);
    if (!c) { printf("ERROR: cut %s missing in %s\n", n, f->GetName()); return nullptr; }
    return c;
}

double M4(const RIn& r) {
    double E = 0, px = 0, py = 0, pz = 0;
    for (int k = 0; k < kNP; ++k) {
        const double p = r.x[k].pc / 1000., th = r.x[k].th * M_PI / 180., ph = r.x[k].ph * M_PI / 180., m = kMass[k] / 1000.;
        px += p * std::sin(th) * std::cos(ph); py += p * std::sin(th) * std::sin(ph); pz += p * std::cos(th);
        E += std::sqrt(p * p + m * m);
    }
    const double m2 = E * E - px * px - py * py - pz * pz;
    return m2 > 0 ? std::sqrt(m2) : 0.;
}

// isData: stored = one segment; sim: stored = two segments
std::vector<Row> Process(const char* file, const char* cutfile, bool isData, Verify& V, const char* tag) {
    std::vector<Row> rows;
    TFile* F = TFile::Open(file);
    TTree* t = (TTree*)F->Get("PipPimEpEm_ID");
    RIn r;
    t->SetBranchStatus("*", 0);
    auto on = [&](const std::string& n, float* a) { t->SetBranchStatus(n.c_str(), 1); t->SetBranchAddress(n.c_str(), a); };
    for (int k = 0; k < kNP; ++k) {
        const std::string s = kPart[k];
        on(s + "_p", &r.x[k].p); on(s + "_" + kCorr[k], &r.x[k].pc); on(s + "_beta", &r.x[k].beta);
        on(s + "_tof", &r.x[k].tof); on(s + "_tof_mom", &r.x[k].tm); on(s + "_tof_new", &r.x[k].tn);
        on(s + "_beta_new", &r.x[k].bn); on(s + "_track_length", &r.x[k].tl); on(s + "_tracklength", &r.x[k].L);
        on(s + "_resolution", &r.x[k].res); on(s + "_theta", &r.x[k].th); on(s + "_phi", &r.x[k].ph);
    }
    on("isBest", &r.isBest); on("runnumber", &r.run); on("event", &r.ev); on("trigbit", &r.trig);
    on("eVertReco_z", &r.vz); on("em_hyp_best_chi2_red", &r.chi2);

    TFile* CF = TFile::Open(cutfile);
    TCutG* cut[kNP];
    for (int k = 0; k < kNP; ++k) cut[k] = GetCut(CF, (std::string(kPart[k]) + (k < 2 ? "_cut_3sig" : "_cut_5sig")).c_str());

    V.hlog[0] = new TH1D(Form("hlog0_%s", tag), "", 140, -10, 0); V.hlog[0]->SetDirectory(0);
    V.hlog[1] = new TH1D(Form("hlog1_%s", tag), "", 140, -10, 0); V.hlog[1]->SetDirectory(0);
    V.hlogbn  = new TH1D(Form("hlogbn_%s", tag), "", 140, -10, 0); V.hlogbn->SetDirectory(0);
    auto lg = [](double d) { return d > 0 ? std::max(-9.95, std::log10(d)) : -9.95; };

    const Long64_t n = t->GetEntries();
    rows.reserve(n);
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        Row o;
        o.run = (long long)r.run; o.ev = (long long)r.ev; o.vz = r.vz; o.isb = r.isBest;
        o.best = (r.isBest == 1);
        const int tb = (int)r.trig;
        o.sel = (r.vz > -160.) && (isData ? (tb == 8192) : ((tb & 8192) != 0));
        double tt[kNP], L[kNP], res[kNP], tmS[kNP], tm1[kNP], tm2[kNP], tmA[kNP];
        for (int k = 0; k < kNP; ++k) {
            const PIn& x = r.x[k];
            tt[k] = x.tof; L[k] = x.L; res[k] = x.res; tmS[k] = x.tm;
            tm1[k] = Tm1(x, kMass[k]); tm2[k] = Tm2(x, kMass[k]);
            o.dtm[k] = tm2[k] - tm1[k];
            o.p[k] = x.p; o.bn0[k] = x.bn;
            const double dOwn = std::fabs(tmS[k] - (isData ? tm1[k] : tm2[k]));
            const double dOth = std::fabs(tmS[k] - (isData ? tm2[k] : tm1[k]));
            if (isData) V.maxd_tm1[k] = std::max(V.maxd_tm1[k], dOwn), V.maxd_tm2[k] = std::max(V.maxd_tm2[k], dOth);
            else        V.maxd_tm2[k] = std::max(V.maxd_tm2[k], dOwn), V.maxd_tm1[k] = std::max(V.maxd_tm1[k], dOth);
            if (dOwn > 1e-4) ++V.nbad_tm;
            V.hlog[0]->Fill(lg(dOwn)); V.hlog[1]->Fill(lg(dOth));
            // alternative = stored + (target - stored formula)
            tmA[k] = isData ? tmS[k] + o.dtm[k] : tmS[k] - o.dtm[k];
        }
        const TC cS = TimeCut(tt, tmS, L, res);   // from the stored tof_mom
        const TC cA = TimeCut(tt, tmA, L, res);   // the other algorithm
        bool rowin = true;
        for (int k = 0; k < kNP; ++k) {
            const double dtn = std::fabs(cS.tn[k] - r.x[k].tn), dbn = std::fabs(cS.bn[k] - r.x[k].bn);
            V.maxd_tn = std::max(V.maxd_tn, dtn); V.maxd_bn = std::max(V.maxd_bn, dbn);
            if (dbn > 1e-5) ++V.nbad_bn;
            V.hlogbn->Fill(lg(dbn));
            o.db[k] = isData ? (cA.bn[k] - cS.bn[k]) : (cS.bn[k] - cA.bn[k]);
            o.in0[k] = cut[k] && cut[k]->IsInside(r.x[k].p, r.x[k].bn);
            const double s = isData ? o.db[k] : -o.db[k];       // change applied to this production
            o.inS[k] = cut[k] && cut[k]->IsInside(r.x[k].p, r.x[k].bn + s);
            o.inM[k] = cut[k] && cut[k]->IsInside(r.x[k].p, r.x[k].bn - s);
            o.inR[k] = o.inRM[k] = o.in0[k];   // filled later (needs the mean shift)
            if (o.in0[k]) ++V.nin0;
            rowin = rowin && o.in0[k];
        }
        if (rowin) ++V.nin_rows;
        const double crel = std::fabs(cS.chi2 - r.chi2) / std::max(1e-6, (double)std::fabs(r.chi2));
        if (std::fabs(cS.chi2 - r.chi2) > 1e-4 && crel > 1e-3) ++V.nbad_chi2;
        if (std::fabs(cS.chi2 - r.chi2) > 1e-4) V.maxd_chi2rel = std::max(V.maxd_chi2rel, crel);
        o.chi2_0 = r.chi2;
        o.chi2_1 = isData ? cS.chi2 : cA.chi2;
        o.chi2_2 = isData ? cA.chi2 : cS.chi2;
        o.m4 = M4(r);
        rows.push_back(o);
        ++V.n;
    }
    V.ntr = 4 * V.n;
    CF->Close(); F->Close();
    return rows;
}

// ---------------------------------------------------------------------------
struct Res { std::vector<double> mu, sig, emu, esig, n; };

void FitBin(TH1D* h, double& mu, double& sig, double& emu, double& esig) {
    mu = h->GetMean(); sig = h->GetRMS(); emu = h->GetMeanError(); esig = h->GetRMSError();
    if (h->GetEntries() < 80) return;
    double m = h->GetBinCenter(h->GetMaximumBin()), s = std::min(h->GetRMS(), 0.03);
    TF1 f("fg", "gaus", m - 2 * s, m + 2 * s);
    for (int it = 0; it < 4; ++it) {
        f.SetParameters(h->GetMaximum(), m, s);
        f.SetRange(m - 2 * s, m + 2 * s);
        h->Fit(&f, "QNR0");
        m = f.GetParameter(1); s = std::fabs(f.GetParameter(2));
    }
    mu = m; sig = s; emu = f.GetParError(1); esig = f.GetParError(2);
}

int Main() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetPalette(kBird);
    gSystem->mkdir(kFig, kTRUE); gSystem->mkdir(kTab, kTRUE);

    Verify VD, VS;
    std::vector<Row> D = Process(kDataFile, kCutExp, true, VD, "d");
    std::vector<Row> S = Process(kSimFile, kCutSim, false, VS, "s");
    printf("data rows %zu, sim rows %zu\n", D.size(), S.size());

    // ------------------------------------------------------------- verification table
    {
        FILE* fo = fopen(Form("%stofmom_verify.txt", kTab), "w");
        for (int s = 0; s < 2; ++s) {
            const Verify& V = s == 0 ? VD : VS;
            fprintf(fo, "# %s file %s, tree PipPimEpEm_ID, %lld rows (%lld particles)\n",
                    s == 0 ? "DATA" : "SIM", s == 0 ? kDataFile : kSimFile, V.n, V.ntr);
            fprintf(fo, "#   stored tof_mom against the %s formula (the production's own):\n", s == 0 ? "ONE-segment" : "TWO-segment");
            for (int k = 0; k < kNP; ++k)
                fprintf(fo, "    %-3s max|stored - one-seg| = %.3e ns   max|stored - two-seg| = %.3e ns\n",
                        kPart[k], V.maxd_tm1[k], V.maxd_tm2[k]);
            fprintf(fo, "    particles with |stored - own formula| > 1e-4 ns: %lld of %lld\n", V.nbad_tm, V.ntr);
            fprintf(fo, "#   tof_new, beta_new, chi2 recomputed with HTimeCut::calcTof from the stored tof, tof_mom, resolution:\n");
            fprintf(fo, "    max|tof_new diff| = %.3e ns, max|beta_new diff| = %.3e, particles with |beta_new diff| > 1e-5: %lld\n",
                    V.maxd_tn, V.maxd_bn, V.nbad_bn);
            fprintf(fo, "    rows with hyp_best_chi2_red differing (abs > 1e-4 and rel > 1e-3): %lld of %lld (max rel %.2e)\n",
                    V.nbad_chi2, V.n, V.maxd_chi2rel);
            fprintf(fo, "#   PID polygons of the production (%s):\n", s == 0 ? "exp betanew_v1" : "sim betanew_v6");
            fprintf(fo, "    particles inside: %lld of %lld; rows with all four inside: %lld of %lld\n\n", V.nin0, V.ntr, V.nin_rows, V.n);
        }
        fclose(fo);
    }

    // ------------------------------------------------------------- resolution per p bin (data, sel, isBest==1)
    Res R[kNP];
    std::vector<TH1D*> hres[kNP];
    for (int k = 0; k < kNP; ++k) {
        const auto& b = Pb(k);
        for (size_t ib = 0; ib + 1 < b.size(); ++ib) {
            TH1D* h = new TH1D(Form("hres_%d_%zu", k, ib), "", 600, -0.15, 0.15); h->SetDirectory(0);
            hres[k].push_back(h);
        }
    }
    for (const Row& o : D) {
        if (!o.sel || !o.best) continue;
        for (int k = 0; k < kNP; ++k) {
            const int ib = PBin(k, o.p[k]); if (ib < 0) continue;
            hres[k][ib]->Fill(o.bn0[k] - BetaP(o.p[k], kMass[k]));
        }
    }
    for (int k = 0; k < kNP; ++k) {
        for (size_t ib = 0; ib < hres[k].size(); ++ib) {
            double mu, sg, emu, esg; FitBin(hres[k][ib], mu, sg, emu, esg);
            R[k].mu.push_back(mu); R[k].sig.push_back(sg); R[k].emu.push_back(emu); R[k].esig.push_back(esg);
            R[k].n.push_back(hres[k][ib]->GetEntries());
        }
    }

    // ------------------------------------------------------------- shift statistics per p bin
    // data: sel rows, isBest==1 (one combination per event); sim the same
    struct Acc { std::vector<double> s0, s1, s2, sa; std::vector<std::vector<double>> v; };   // sa: sum of delta tof_mom
    auto shiftStats = [&](const std::vector<Row>& X, int k, bool onlyBest) {
        const auto& b = Pb(k); const size_t nb = b.size() - 1;
        Acc a; a.s0.assign(nb, 0); a.s1.assign(nb, 0); a.s2.assign(nb, 0); a.sa.assign(nb, 0); a.v.assign(nb, {});
        for (const Row& o : X) {
            if (!o.sel || (onlyBest && !o.best)) continue;
            const int ib = PBin(k, o.p[k]); if (ib < 0) continue;
            a.s0[ib] += 1; a.s1[ib] += o.db[k]; a.s2[ib] += o.db[k] * o.db[k]; a.sa[ib] += o.dtm[k]; a.v[ib].push_back(o.db[k]);
        }
        return a;
    };
    Acc AD[kNP], AS[kNP], ADall[kNP];
    for (int k = 0; k < kNP; ++k) { AD[k] = shiftStats(D, k, true); AS[k] = shiftStats(S, k, true); ADall[k] = shiftStats(D, k, false); }
    auto quant = [](std::vector<double> v, double q) {
        if (v.empty()) return 0.; std::sort(v.begin(), v.end());
        return v[std::min(v.size() - 1, size_t(q * v.size()))];
    };
    // mean shift per bin from all sel rows (for the re-centred variant)
    std::vector<double> meanShift[kNP];
    for (int k = 0; k < kNP; ++k)
        for (size_t ib = 0; ib < ADall[k].s0.size(); ++ib)
            meanShift[k].push_back(ADall[k].s0[ib] > 0 ? ADall[k].s1[ib] / ADall[k].s0[ib] : 0.);
    // polygon status with the mean shift of the bin removed (data only)
    {
        TFile* CF = TFile::Open(kCutExp);
        TCutG* cut[kNP];
        for (int k = 0; k < kNP; ++k) cut[k] = GetCut(CF, (std::string(kPart[k]) + (k < 2 ? "_cut_3sig" : "_cut_5sig")).c_str());
        for (Row& o : D)
            for (int k = 0; k < kNP; ++k) {
                const int ib = PBin(k, o.p[k]);
                const double ms = ib >= 0 ? meanShift[k][ib] : 0.;
                o.inR[k]  = cut[k]->IsInside(o.p[k], o.bn0[k] + o.db[k] - ms);
                o.inRM[k] = cut[k]->IsInside(o.p[k], o.bn0[k] - o.db[k] + ms);
            }
        CF->Close();
    }

    {
        FILE* fo = fopen(Form("%stofmom_shift_vs_p.txt", kTab), "w");
        fprintf(fo, "# delta beta_new = beta_new(two-segment tof_mom) - beta_new(one-segment), per particle and p bin\n");
        fprintf(fo, "# DATA %s (sel: eVertReco_z > -160, trigbit == 8192, isBest == 1); SIM %s (sel: vz > -160, trigbit & 8192, isBest == 1)\n", kDataFile, kSimFile);
        fprintf(fo, "# sigma = measured beta_new resolution of the data: gaussian core (+-2 sigma, iterated) of beta_new - beta(p, hypothesis mass)\n");
        fprintf(fo, "# columns: p range [MeV/c] | data N | data mean dtof_mom [ns] | data mean dbeta | rms | 95%% q of dbeta | data mu | sigma (+-err) | mean dbeta / sigma | 95%% q / sigma | sim N | sim mean dbeta | sim rms\n");
        for (int k = 0; k < kNP; ++k) {
            fprintf(fo, "\n## %s\n", kPart[k]);
            const auto& b = Pb(k);
            for (size_t ib = 0; ib + 1 < b.size(); ++ib) {
                const double N = AD[k].s0[ib], m = N > 0 ? AD[k].s1[ib] / N : 0, rms = N > 1 ? std::sqrt(std::max(0., AD[k].s2[ib] / N - m * m)) : 0;
                const double q95 = quant(AD[k].v[ib], 0.95);
                const double NS = AS[k].s0[ib], mS = NS > 0 ? AS[k].s1[ib] / NS : 0, rS = NS > 1 ? std::sqrt(std::max(0., AS[k].s2[ib] / NS - mS * mS)) : 0;
                const double sg = R[k].sig[ib];
                fprintf(fo, "%5.0f-%5.0f | %7.0f | %+.4f | %+.5f | %.5f | %+.5f | %+.4f | %.4f +- %.4f | %+.3f | %+.3f | %7.0f | %+.5f | %.5f\n",
                        b[ib], b[ib + 1], N, N > 0 ? AD[k].sa[ib] / N : 0., m, rms, q95, R[k].mu[ib], sg, R[k].esig[ib], sg > 0 ? m / sg : 0., sg > 0 ? q95 / sg : 0., NS, mS, rS);
            }
        }
        fclose(fo);
    }

    // ------------------------------------------------------------- crossing counts (data, sel, all rows = isBest >= 0)
    struct Cnt { double n = 0, loss = 0, gain = 0, lossR = 0, gainR = 0, lossN = 0, gainN = 0, lossNR = 0, gainNR = 0; };
    Cnt cp[kNP], crow, crowBest;
    auto nsig = [&](int k, double p, double bn, double& z) {
        const int ib = PBin(k, p); if (ib < 0 || R[k].sig[ib] <= 0) return false;
        z = (bn - BetaP(p, kMass[k]) - R[k].mu[ib]) / R[k].sig[ib]; return true;
    };
    std::vector<double> lossVsP[kNP], gainVsP[kNP], allVsP[kNP], lossRVsP[kNP], gainRVsP[kNP];
    for (int k = 0; k < kNP; ++k) { const size_t nb = Pb(k).size() - 1; lossVsP[k].assign(nb, 0); gainVsP[k].assign(nb, 0); allVsP[k].assign(nb, 0);
                                    lossRVsP[k].assign(nb, 0); gainRVsP[k].assign(nb, 0); }
    TH1D* hm4[5];
    const char* hm4n[5] = {"all", "loss", "gain", "lossR", "gainR"};
    for (int i = 0; i < 5; ++i) { hm4[i] = new TH1D(Form("hm4_%s", hm4n[i]), "", 40, 0.30, 1.10); hm4[i]->Sumw2(); hm4[i]->SetDirectory(0); }
    for (const Row& o : D) {
        if (!o.sel) continue;
        bool rl = false, rg = false, rlR = false, rgR = false, rlN = false, rgN = false, rlNR = false, rgNR = false;
        bool rowin = true;
        for (int k = 0; k < kNP; ++k) rowin = rowin && o.in0[k];
        if (!rowin) continue;   // (every row is inside, see verification)
        for (int k = 0; k < kNP; ++k) {
            Cnt& c = cp[k]; c.n += 1;
            const bool l = o.in0[k] && !o.inS[k], g = o.in0[k] && !o.inM[k];
            const bool lR = o.in0[k] && !o.inR[k], gR = o.in0[k] && !o.inRM[k];
            // mirror of the re-centred shift
            const int ib = PBin(k, o.p[k]);
            const double ms = ib >= 0 ? meanShift[k][ib] : 0.;
            double z0 = 0, z1 = 0, zm = 0, zr = 0, zrm = 0;
            bool lN = false, gN = false, lNR = false, gNR = false;
            if (nsig(k, o.p[k], o.bn0[k], z0) && std::fabs(z0) < kNsig[k]) {
                nsig(k, o.p[k], o.bn0[k] + o.db[k], z1); nsig(k, o.p[k], o.bn0[k] - o.db[k], zm);
                nsig(k, o.p[k], o.bn0[k] + o.db[k] - ms, zr); nsig(k, o.p[k], o.bn0[k] - o.db[k] + ms, zrm);
                lN = std::fabs(z1) >= kNsig[k]; gN = std::fabs(zm) >= kNsig[k];
                lNR = std::fabs(zr) >= kNsig[k]; gNR = std::fabs(zrm) >= kNsig[k];
            }
            c.loss += l; c.gain += g; c.lossR += lR; c.gainR += gR; c.lossN += lN; c.gainN += gN; c.lossNR += lNR; c.gainNR += gNR;
            rl |= l; rg |= g; rlR |= lR; rgR |= gR; rlN |= lN; rgN |= gN; rlNR |= lNR; rgNR |= gNR;
            if (ib >= 0) { allVsP[k][ib] += 1; lossVsP[k][ib] += l; gainVsP[k][ib] += g; lossRVsP[k][ib] += lR; gainRVsP[k][ib] += gR; }
        }
        // mirror gain for the re-centred variant (polygon): counted per row below
        for (Cnt* c : {&crow, &crowBest}) {
            if (c == &crowBest && !o.best) continue;
            c->n += 1; c->loss += rl; c->gain += rg; c->lossR += rlR; c->gainR += rgR;
            c->lossN += rlN; c->gainN += rgN; c->lossNR += rlNR; c->gainNR += rgNR;
        }
        hm4[0]->Fill(o.m4); if (rl) hm4[1]->Fill(o.m4); if (rg) hm4[2]->Fill(o.m4); if (rlR) hm4[3]->Fill(o.m4); if (rgR) hm4[4]->Fill(o.m4);
    }
    // ------------------------------------------------------------- isBest (data, all rows; groups = adjacent equal (run, event))
    long long ngrp = 0, ngrpMulti = 0, nStoredOk = 0, nChange = 0, nChangeSel = 0, ngrpMultiSel = 0;
    double chi2Old = 0, chi2New = 0; long long nchi = 0;
    {
        size_t i0 = 0;
        while (i0 < D.size()) {
            size_t i1 = i0 + 1;
            while (i1 < D.size() && D[i1].run == D[i0].run && D[i1].ev == D[i0].ev) ++i1;
            ++ngrp;
            size_t a0 = i0, a1 = i0, a2 = i0; int nb1 = 0;
            for (size_t i = i0; i < i1; ++i) {
                if (D[i].chi2_0 < D[a0].chi2_0) a0 = i;
                if (D[i].chi2_1 < D[a1].chi2_1) a1 = i;
                if (D[i].chi2_2 < D[a2].chi2_2) a2 = i;
                if (D[i].best) ++nb1;
                chi2Old += D[i].chi2_1; chi2New += D[i].chi2_2; ++nchi;
            }
            if (nb1 == 1 && D[a0].best) ++nStoredOk;
            if (i1 - i0 > 1) {
                ++ngrpMulti; if (D[i0].sel) ++ngrpMultiSel;
                if (a2 != a1) { ++nChange; if (D[i0].sel) ++nChangeSel; }
            }
            D[a2].newBest = true;
            i0 = i1;
        }
    }

    {
        FILE* fo = fopen(Form("%stofmom_crossing.txt", kTab), "w");
        fprintf(fo, "# DATA %s, PipPimEpEm_ID, sel: eVertReco_z > -160, trigbit == 8192, all rows (isBest >= 0)\n", kDataFile);
        fprintf(fo, "# change: tof_mom one segment -> two segments (HTimeCut of the simulation), tof_new/beta_new recomputed with HTimeCut::calcTof\n");
        fprintf(fo, "# loss  = inside now, outside after the change; gain = mirror estimate of the tracks entering from outside\n");
        fprintf(fo, "#         (= tracks inside within the shift of the opposite edge; the outside tracks are not in the _ID file)\n");
        fprintf(fo, "# POLY  = data PID polygons pp45_pid_cuts_exp_betanew_v1.root (pions 3sig, leptons 5sig) at the data p\n");
        fprintf(fo, "# NSIG  = |beta_new - beta(p) - mu(p)| < k sigma(p) with the measured mu, sigma (k = 3 pions, 5 leptons)\n");
        fprintf(fo, "# RECENT= the same after removing the mean shift of the p bin (polygons refitted after a reprocessing)\n");
        fprintf(fo, "# fractions in %% of the particles (rows); binomial errors (tracks repeat in several rows: underestimated)\n\n");
        auto pr = [&](const char* name, const Cnt& c) {
            auto f = [&](double x) { return 100. * x / std::max(1., c.n); };
            auto e = [&](double x) { const double q = x / std::max(1., c.n); return 100. * std::sqrt(q * (1 - q) / std::max(1., c.n)); };
            fprintf(fo, "%-10s N=%8.0f | POLY loss %.3f+-%.3f gain %.3f+-%.3f net %+.3f | NSIG loss %.3f gain %.3f net %+.3f | POLY RECENT loss %.3f gain %.3f net %+.3f | NSIG RECENT loss %.3f gain %.3f net %+.3f\n",
                    name, c.n, f(c.loss), e(c.loss), f(c.gain), e(c.gain), f(c.gain) - f(c.loss), f(c.lossN), f(c.gainN), f(c.gainN) - f(c.lossN),
                    f(c.lossR), f(c.gainR), f(c.gainR) - f(c.lossR), f(c.lossNR), f(c.gainNR), f(c.gainNR) - f(c.lossNR));
        };
        for (int k = 0; k < kNP; ++k) pr(kPart[k], cp[k]);
        pr("rows", crow);
        pr("rows best", crowBest);
        fprintf(fo, "\n# loss / gain per p bin (POLY and POLY RECENT), %% of the particles of the bin\n");
        for (int k = 0; k < kNP; ++k) {
            fprintf(fo, "## %s\n", kPart[k]);
            const auto& b = Pb(k);
            for (size_t ib = 0; ib + 1 < b.size(); ++ib)
                fprintf(fo, "%5.0f-%5.0f N=%7.0f loss %.3f gain %.3f | RECENT loss %.3f gain %.3f\n", b[ib], b[ib + 1], allVsP[k][ib],
                        allVsP[k][ib] > 0 ? 100. * lossVsP[k][ib] / allVsP[k][ib] : 0., allVsP[k][ib] > 0 ? 100. * gainVsP[k][ib] / allVsP[k][ib] : 0.,
                        allVsP[k][ib] > 0 ? 100. * lossRVsP[k][ib] / allVsP[k][ib] : 0., allVsP[k][ib] > 0 ? 100. * gainRVsP[k][ib] / allVsP[k][ib] : 0.);
        }
        fprintf(fo, "\n# isBest (HTimeCut alpha = 0: lowest hyp_best_chi2_red among the PipPimEpEm rows of an event), all rows of the file\n");
        fprintf(fo, "events %lld, with >= 2 rows %lld (sel %lld); stored isBest==1 is the stored chi2 minimum in %lld events (%.3f %%)\n",
                ngrp, ngrpMulti, ngrpMultiSel, nStoredOk, 100. * nStoredOk / std::max(1LL, ngrp));
        fprintf(fo, "best row changes with the two-segment chi2: %lld of %lld multi-row events (%.2f %%); sel %lld of %lld (%.2f %%)\n",
                nChange, ngrpMulti, 100. * nChange / std::max(1LL, ngrpMulti), nChangeSel, ngrpMultiSel, 100. * nChangeSel / std::max(1LL, ngrpMultiSel));
        fprintf(fo, "mean hyp_best_chi2_red one segment %.4f, two segments %.4f\n", chi2Old / std::max(1LL, nchi), chi2New / std::max(1LL, nchi));
        // M4 regions
        fprintf(fo, "\n# rows (sel) by M4 region (M4 from p_corr), %% of N; POLY = data polygons kept, RECENT = mean shift removed;\n");
        fprintf(fo, "# isBest==1: rows that are best with the one-segment chi2 (stored) and with the two-segment chi2 (PID unchanged)\n");
        fprintf(fo, "# region       range         | N(isBest>=0) | POLY loss gain net (+-) | RECENT loss gain net (+-) | N(isBest==1) old new change%%\n");
        const double reg[6][2] = {{0.30, 0.455}, {0.455, 0.5225}, {0.5336, 0.5542}, {0.58, 0.70}, {0.70, 1.10}, {0.30, 1.10}};
        const char* regn[6] = {"below bump", "bump", "peak", "above peak", "high mass", "all"};
        for (int ir = 0; ir < 6; ++ir) {
            double N = 0, l = 0, g = 0, lR = 0, gR = 0, nb0 = 0, nb1 = 0;
            for (const Row& o : D) {
                if (!o.sel || o.m4 < reg[ir][0] || o.m4 >= reg[ir][1]) continue;
                bool rl = false, rg = false, rlR = false, rgR = false;
                for (int k = 0; k < kNP; ++k) { rl |= !o.inS[k]; rg |= !o.inM[k]; rlR |= !o.inR[k]; rgR |= !o.inRM[k]; }
                N += 1; l += rl; g += rg; lR += rlR; gR += rgR; nb0 += o.best; nb1 += o.newBest;
            }
            const double n1 = std::max(1., N);
            fprintf(fo, "%-11s %.4f-%.4f | %7.0f | %.3f %.3f %+.3f (%.3f) | %.3f %.3f %+.3f (%.3f) | %6.0f %6.0f %+.2f\n", regn[ir], reg[ir][0], reg[ir][1], N,
                    100. * l / n1, 100. * g / n1, 100. * (g - l) / n1, 100. * std::sqrt(l + g) / n1,
                    100. * lR / n1, 100. * gR / n1, 100. * (gR - lR) / n1, 100. * std::sqrt(lR + gR) / n1,
                    nb0, nb1, nb0 > 0 ? 100. * (nb1 - nb0) / nb0 : 0.);
        }
        {   // simulation chi2 means
            double c1 = 0, c2 = 0; long long nn = 0;
            for (const Row& o : S) { c1 += o.chi2_1; c2 += o.chi2_2; ++nn; }
            fprintf(fo, "\nSIM mean hyp_best_chi2_red: one segment %.4f, two segments (stored) %.4f\n", c1 / std::max(1LL, nn), c2 / std::max(1LL, nn));
        }
        fclose(fo);
    }

    // ------------------------------------------------------------- simulation: reverse change on its own polygons
    {
        Cnt sp[kNP], srow;
        for (const Row& o : S) {
            if (!o.sel) continue;
            bool rowin = true; for (int k = 0; k < kNP; ++k) rowin = rowin && o.in0[k];
            if (!rowin) continue;
            bool rl = false, rg = false;
            for (int k = 0; k < kNP; ++k) {
                const bool l = !o.inS[k], g = !o.inM[k];
                sp[k].n += 1; sp[k].loss += l; sp[k].gain += g; rl |= l; rg |= g;
            }
            srow.n += 1; srow.loss += rl; srow.gain += rg;
        }
        FILE* fo = fopen(Form("%stofmom_crossing.txt", kTab), "a");
        fprintf(fo, "\n# SIM %s (unweighted rows; sel vz > -160, trigbit & 8192, isBest >= 0): change two -> one segment on the sim polygons v6\n", kSimFile);
        for (int k = 0; k < kNP; ++k)
            fprintf(fo, "%-4s N=%8.0f loss %.3f gain %.3f net %+.3f %%\n", kPart[k], sp[k].n, 100. * sp[k].loss / std::max(1., sp[k].n),
                    100. * sp[k].gain / std::max(1., sp[k].n), 100. * (sp[k].gain - sp[k].loss) / std::max(1., sp[k].n));
        fprintf(fo, "rows N=%8.0f loss %.3f gain %.3f net %+.3f %%\n", srow.n, 100. * srow.loss / std::max(1., srow.n),
                100. * srow.gain / std::max(1., srow.n), 100. * (srow.gain - srow.loss) / std::max(1., srow.n));
        fclose(fo);
    }

    // ============================================================= figures
    auto txt = [](double x, double y, const char* s, double size, int col = 1, int al = 11) {
        TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(al); t.DrawLatex(x, y, s);
    };
    auto ax = [](TH1* h, const char* xt, const char* yt, double sc = 1.) {
        h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
        h->GetXaxis()->SetTitleSize(0.055 * sc); h->GetYaxis()->SetTitleSize(0.055 * sc);
        h->GetXaxis()->SetLabelSize(0.048 * sc); h->GetYaxis()->SetLabelSize(0.048 * sc);
        h->GetXaxis()->SetTitleOffset(1.05); h->GetYaxis()->SetTitleOffset(1.3 / sc);
        h->GetYaxis()->SetMaxDigits(3);
    };

    // (1) verification
    {
        TCanvas c("cv", "", 1500, 600); c.Divide(2, 1);
        for (int s = 0; s < 2; ++s) {
            const Verify& V = s == 0 ? VD : VS;
            c.cd(s + 1); gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.12);
            TH1D* a = V.hlog[0]; TH1D* b = V.hlog[1]; TH1D* bn = V.hlogbn;
            a->SetLineColor(kBlue + 1); a->SetLineWidth(3); b->SetLineColor(kRed + 1); b->SetLineWidth(3); b->SetLineStyle(2);
            bn->SetLineColor(kBlack); bn->SetLineWidth(2);
            a->SetMaximum(300 * std::max({a->GetMaximum(), b->GetMaximum(), bn->GetMaximum()})); a->SetMinimum(0.5);
            ax(a, "log_{10}|stored - recomputed|  (-10 = identical)", "particles");
            a->Draw("hist"); b->Draw("hist same"); bn->Draw("hist same");
            TLegend l(0.15, 0.68, 0.88, 0.87); l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextSize(0.033);
            l.AddEntry(a, s == 0 ? "t_{mom} [ns] vs one-segment formula (data PAT): reproduced" : "t_{mom} [ns] vs two-segment formula (sim PAT): reproduced", "l");
            l.AddEntry(b, s == 0 ? "t_{mom} [ns] vs two-segment formula (not used in data)" : "t_{mom} [ns] vs one-segment formula (not used in sim)", "l");
            l.AddEntry(bn, "#beta_{new} vs HTimeCut::calcTof from stored t, t_{mom}", "l");
            l.DrawClone();
            txt(0.13, 0.91, Form("%s, PipPimEpEm_ID: %lld rows x 4 particles", s == 0 ? "DATA lep_034" : "SIM ppee_std_all", V.n), 0.040);
        }
        c.SaveAs(Form("%stofmom_verify.png", kFig));
    }

    // (2) delta beta_new vs p, data, per particle with the measured sigma
    std::vector<TGraphErrors*> gSig(kNP), gMean(kNP), gMeanS(kNP), gRatio(kNP), gQ95(kNP);
    for (int k = 0; k < kNP; ++k) {
        const auto& b = Pb(k);
        gSig[k] = new TGraphErrors(); gMean[k] = new TGraphErrors(); gMeanS[k] = new TGraphErrors(); gRatio[k] = new TGraphErrors(); gQ95[k] = new TGraphErrors();
        for (size_t ib = 0; ib + 1 < b.size(); ++ib) {
            const double pc = 0.5 * (b[ib] + b[ib + 1]), ex = 0.5 * (b[ib + 1] - b[ib]);
            const double N = AD[k].s0[ib]; if (N < 200) continue;
            const double m = AD[k].s1[ib] / N, rms = std::sqrt(std::max(0., AD[k].s2[ib] / N - m * m));
            const double sg = R[k].sig[ib];
            int i = gSig[k]->GetN();
            gSig[k]->SetPoint(i, pc, sg); gSig[k]->SetPointError(i, ex, R[k].esig[ib]);
            gMean[k]->SetPoint(i, pc, m); gMean[k]->SetPointError(i, ex, rms / std::sqrt(N));
            gRatio[k]->SetPoint(i, pc, m / sg); gRatio[k]->SetPointError(i, ex, std::hypot(rms / std::sqrt(N) / sg, m * R[k].esig[ib] / (sg * sg)));
            gQ95[k]->SetPoint(i, pc, quant(AD[k].v[ib], 0.95) / sg); gQ95[k]->SetPointError(i, ex, 0);
            const double NS = AS[k].s0[ib];
            if (NS >= 20) { const double mS = AS[k].s1[ib] / NS, rS = std::sqrt(std::max(0., AS[k].s2[ib] / NS - mS * mS));
                int j = gMeanS[k]->GetN(); gMeanS[k]->SetPoint(j, pc, mS); gMeanS[k]->SetPointError(j, ex, rS / std::sqrt(NS)); }
        }
    }
    {
        TCanvas c("cd", "", 1500, 1100); c.Divide(2, 2);
        TH2D* h2[kNP];
        for (int k = 0; k < kNP; ++k) {
            const bool pi = k < 2;
            h2[k] = new TH2D(Form("h2db_%d", k), "", pi ? 62 : 48, pi ? 60 : 0, pi ? 1300 : 1200, 160, pi ? -0.002 : -0.003, pi ? 0.030 : 0.013);
            h2[k]->SetDirectory(0);
        }
        for (const Row& o : D) { if (!o.sel || !o.best) continue; for (int k = 0; k < kNP; ++k) h2[k]->Fill(o.p[k], o.db[k]); }
        for (int k = 0; k < kNP; ++k) {
            c.cd(k + 1); gPad->SetLogz(); gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.12); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.14);
            ax(h2[k], Form("p_{%s} (track fit) [MeV/c]", kLabel[k]), Form("#Delta#beta_{new}(%s) = two-seg. - one-seg.", kLabel[k]));
            h2[k]->Draw("colz");
            gSig[k]->SetLineColor(kRed + 1); gSig[k]->SetLineWidth(3); gSig[k]->SetMarkerColor(kRed + 1); gSig[k]->SetMarkerStyle(20);
            gSig[k]->Draw("LP same");
            gMean[k]->SetMarkerStyle(20); gMean[k]->SetMarkerSize(1.2); gMean[k]->SetMarkerColor(kBlack); gMean[k]->SetLineColor(kBlack);
            gMean[k]->Draw("P same");
            TLine z(h2[k]->GetXaxis()->GetXmin(), 0, h2[k]->GetXaxis()->GetXmax(), 0); z.SetLineStyle(2); z.DrawClone();
            TLegend l(0.42, 0.72, 0.87, 0.89); l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextSize(0.042);
            l.AddEntry(gMean[k], "mean #Delta#beta_{new} per p bin", "p");
            l.AddEntry(gSig[k], "measured #sigma(#beta_{new}) of the data", "lp");
            l.DrawClone();
            txt(0.14, 0.94, Form("DATA lep_034, PipPimEpEm_ID, vz > -160, PT3, isBest = 1: %s", kLabel[k]), 0.042);
        }
        c.SaveAs(Form("%stofmom_dbeta_vs_p.png", kFig));
        c.SaveAs(Form("%stofmom_dbeta_vs_p.pdf", kFig));
    }

    // (3) shift in units of sigma
    {
        TCanvas c("cr", "", 1500, 650); c.Divide(2, 1);
        for (int half = 0; half < 2; ++half) {
            c.cd(half + 1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.1); gPad->SetLogx(half == 0);
            TH1D fr(Form("fr%d", half), "", 10, half == 0 ? 60 : 0, half == 0 ? 2500 : 1900);
            fr.SetMinimum(half == 0 ? -0.05 : -0.05); fr.SetMaximum(half == 0 ? 1.2 : 0.8);
            ax(&fr, "p (track fit) [MeV/c]", "shift / #sigma");
            fr.GetXaxis()->SetNdivisions(505);
            fr.DrawCopy();
            TLegend l(0.45, 0.62, 0.88, 0.88); l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextSize(0.04);
            const int col[kNP] = {kRed + 1, kBlue + 1, kMagenta + 1, kGreen + 2};
            for (int k = half * 2; k < half * 2 + 2; ++k) {
                gRatio[k]->SetMarkerStyle(20 + k % 2); gRatio[k]->SetMarkerColor(col[k]); gRatio[k]->SetLineColor(col[k]); gRatio[k]->SetMarkerSize(1.3);
                gRatio[k]->Draw("P same");
                gQ95[k]->SetLineColor(col[k]); gQ95[k]->SetLineStyle(2); gQ95[k]->SetLineWidth(2);
                gQ95[k]->Draw("L same");
                l.AddEntry(gRatio[k], Form("%s: mean shift", kLabel[k]), "p");
                l.AddEntry(gQ95[k], Form("%s: 95%% quantile of the shift", kLabel[k]), "l");
            }
            l.DrawClone();
            TLine z(half == 0 ? 60 : 0, 0, half == 0 ? 2500 : 1900, 0); z.SetLineStyle(3); z.DrawClone();
            txt(0.13, 0.92, half == 0 ? "DATA pions: #Delta#beta_{new} / measured #sigma(#beta_{new}), isBest = 1"
                                      : "DATA leptons: the same (e^{-} = time reference)", 0.036);
        }
        c.SaveAs(Form("%stofmom_shift_sigma.png", kFig));
        c.SaveAs(Form("%stofmom_shift_sigma.pdf", kFig));
    }

    // (4) data vs simulation mean shift
    {
        TCanvas c("cs", "", 1500, 650); c.Divide(2, 1);
        for (int half = 0; half < 2; ++half) {
            c.cd(half + 1); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.14);
            TH1D fr(Form("fs%d", half), "", 10, half == 0 ? 60 : 0, half == 0 ? 1300 : 1200);
            fr.SetMinimum(half == 0 ? -0.0005 : -0.0005); fr.SetMaximum(half == 0 ? 0.012 : 0.004);
            ax(&fr, "p (track fit) [MeV/c]", "mean #Delta#beta_{new} (two-seg. - one-seg.)");
            fr.DrawCopy();
            TLegend l(0.45, 0.62, 0.88, 0.88); l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextSize(0.04);
            const int col[kNP] = {kRed + 1, kBlue + 1, kMagenta + 1, kGreen + 2};
            for (int k = half * 2; k < half * 2 + 2; ++k) {
                TGraphErrors* gd = (TGraphErrors*)gMean[k]->Clone(); gd->SetMarkerColor(col[k]); gd->SetLineColor(col[k]); gd->SetMarkerStyle(20);
                gd->Draw("P same");
                gMeanS[k]->SetLineColor(col[k]); gMeanS[k]->SetLineWidth(3); gMeanS[k]->SetLineStyle(2); gMeanS[k]->SetMarkerSize(0);
                gMeanS[k]->Draw("L same");
                l.AddEntry(gd, Form("%s data (two-seg. would add)", kLabel[k]), "p");
                l.AddEntry(gMeanS[k], Form("%s simulation (two-seg. added)", kLabel[k]), "l");
            }
            l.DrawClone();
            txt(0.15, 0.94, "Mean change of #beta_{new} per p bin; vz > -160, PT3, isBest = 1", 0.040);
        }
        c.SaveAs(Form("%stofmom_data_vs_sim_shift.png", kFig));
    }

    // (5) crossing fraction vs p (polygon)
    {
        TCanvas c("cx", "", 1500, 1100); c.Divide(2, 2);
        for (int k = 0; k < kNP; ++k) {
            c.cd(k + 1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.1);
            const auto& b = Pb(k);
            TH1D* hh[3]; const char* nm[3] = {"hl", "hg", "hlr"};
            std::vector<double>* src[3] = {lossVsP, gainVsP, lossRVsP};
            for (int j = 0; j < 3; ++j) {
                hh[j] = new TH1D(Form("%s%d", nm[j], k), "", b.size() - 1, b.data()); hh[j]->SetDirectory(0);
                for (size_t ib = 0; ib + 1 < b.size(); ++ib) {
                    const double N = allVsP[k][ib]; if (N < 20) continue;
                    const double f = src[j][k][ib] / N;
                    hh[j]->SetBinContent(ib + 1, 100 * f); hh[j]->SetBinError(ib + 1, 100 * std::sqrt(std::max(f, 1. / N) * (1 - f) / N));
                }
            }
            ax(hh[0], Form("p_{%s} (track fit) [MeV/c]", kLabel[k]), "tracks crossing the polygon edge [%]");
            hh[0]->GetXaxis()->SetNdivisions(505);
            const int col[3] = {kRed + 1, kBlue + 1, kGreen + 2}, mk[3] = {20, 24, 21};
            for (int j = 0; j < 3; ++j) { hh[j]->SetLineColor(col[j]); hh[j]->SetMarkerColor(col[j]); hh[j]->SetMarkerStyle(mk[j]); hh[j]->SetLineWidth(2); }
            if (k < 2) { gPad->SetLogx(); gPad->SetLogy(); hh[0]->SetMinimum(0.005); hh[0]->SetMaximum(300); }
            else { hh[0]->SetMinimum(0); hh[0]->SetMaximum(0.5); }
            hh[0]->Draw("E1"); hh[1]->Draw("E1 same"); hh[2]->Draw("E1 same");
            TLegend l(0.36, 0.70, 0.88, 0.88); l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextSize(0.038);
            l.AddEntry(hh[0], "leave the band, data polygons kept", "p");
            l.AddEntry(hh[1], "enter the band (mirror estimate)", "p");
            l.AddEntry(hh[2], "leave, mean shift removed (polygons refitted)", "p");
            l.DrawClone();
            if (k == 2) txt(0.30, 0.45, "e^{+}: #Delta#beta_{new} < 10^{-6}, no track crosses", 0.045);
            txt(0.14, 0.92, Form("DATA: %s, data polygon %s; vz > -160, PT3, isBest #geq 0", kLabel[k], k < 2 ? "3#sigma" : "5#sigma"), 0.040);
        }
        c.SaveAs(Form("%stofmom_crossing_vs_p.png", kFig));
        c.SaveAs(Form("%stofmom_crossing_vs_p.pdf", kFig));
    }

    // (6) M4: rows lost and gained
    {
        TCanvas c("cm", "", 1400, 650); c.Divide(2, 1);
        c.cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.1);
        ax(hm4[0], "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (p_{corr}) [GeV/c^{2}]", "rows / 20 MeV/c^{2}");
        hm4[0]->SetMarkerStyle(20); hm4[0]->SetMarkerColor(kBlack); hm4[0]->SetLineColor(kBlack);
        hm4[1]->SetMarkerStyle(21); hm4[1]->SetMarkerColor(kRed + 1); hm4[1]->SetLineColor(kRed + 1);
        hm4[2]->SetMarkerStyle(25); hm4[2]->SetMarkerColor(kBlue + 1); hm4[2]->SetLineColor(kBlue + 1);
        hm4[0]->SetMinimum(0.5); hm4[0]->SetMaximum(300 * hm4[0]->GetMaximum());
        hm4[0]->Draw("E1"); hm4[1]->Draw("E1 same"); hm4[2]->Draw("E1 same");
        TLegend l(0.17, 0.70, 0.88, 0.88); l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextSize(0.04);
        l.AddEntry(hm4[0], "all rows (data, one segment)", "p");
        l.AddEntry(hm4[1], "rows that leave (two segments, data polygons)", "p");
        l.AddEntry(hm4[2], "rows that enter (mirror estimate)", "p");
        l.DrawClone();
        txt(0.14, 0.92, "DATA lep_034; vz > -160, PT3, isBest #geq 0", 0.042);
        c.cd(2); gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.1);
        TH1D* hr = (TH1D*)hm4[0]->Clone("hm4_ratio"); hr->Reset();
        TH1D* hrR = (TH1D*)hm4[0]->Clone("hm4_ratioR"); hrR->Reset();
        for (int i = 1; i <= hr->GetNbinsX(); ++i) {
            const double N = hm4[0]->GetBinContent(i); if (N <= 0) continue;
            const double lo = hm4[1]->GetBinContent(i), ga = hm4[2]->GetBinContent(i);
            const double loR = hm4[3]->GetBinContent(i), gaR = hm4[4]->GetBinContent(i);
            hr->SetBinContent(i, (N - lo + ga) / N); hr->SetBinError(i, std::sqrt(lo + ga) / N);
            hrR->SetBinContent(i, (N - loR + gaR) / N); hrR->SetBinError(i, std::sqrt(loR + gaR) / N);
        }
        ax(hr, "M_{#pi^{+}#pi^{-}e^{+}e^{-}} (p_{corr}) [GeV/c^{2}]", "rows after / rows before");
        hr->SetMinimum(0.97); hr->SetMaximum(1.02); hr->GetYaxis()->SetTitleOffset(1.6); hr->SetMarkerStyle(20); hr->SetMarkerColor(kRed + 1); hr->SetLineColor(kRed + 1);
        hrR->SetMarkerStyle(21); hrR->SetMarkerColor(kGreen + 2); hrR->SetLineColor(kGreen + 2);
        hr->Draw("E1"); hrR->Draw("E1 same");
        TLine one(0.30, 1, 1.10, 1); one.SetLineStyle(2); one.DrawClone();
        TLegend l2(0.20, 0.16, 0.88, 0.30); l2.SetBorderSize(0); l2.SetFillStyle(0); l2.SetTextSize(0.038);
        l2.AddEntry(hr, "two segments, data polygons kept", "p");
        l2.AddEntry(hrR, "two segments, polygons refitted (mean shift removed)", "p");
        l2.DrawClone();
        txt(0.14, 0.92, "(N - leave + enter) / N per 20 MeV/c^{2}; linear zoom", 0.042);
        c.SaveAs(Form("%stofmom_m4.png", kFig));
        c.SaveAs(Form("%stofmom_m4.pdf", kFig));
    }

    // (7) example slice: pi+ and e- beta_new residual in a low-p bin, before/after, with the polygon edges
    {
        TFile* CF = TFile::Open(kCutExp);
        TCanvas c("ce", "", 1500, 650); c.Divide(2, 1);
        const int kk[2] = {0, 3}; const double plo[2] = {100, 100}, phi[2] = {130, 200};
        for (int s = 0; s < 2; ++s) {
            const int k = kk[s];
            TCutG* cut = GetCut(CF, (std::string(kPart[k]) + (k < 2 ? "_cut_3sig" : "_cut_5sig")).c_str());
            c.cd(s + 1); gPad->SetLogy(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.1);
            const double pm = 0.5 * (plo[s] + phi[s]);
            const double b0 = BetaP(pm, kMass[k]);
            TH1D* h0 = new TH1D(Form("hs0_%d", s), "", 200, -0.10, 0.10); TH1D* h1 = new TH1D(Form("hs1_%d", s), "", 200, -0.10, 0.10);
            h0->SetDirectory(0); h1->SetDirectory(0);
            for (const Row& o : D) {
                if (!o.sel || o.p[k] < plo[s] || o.p[k] >= phi[s]) continue;
                const double bp = BetaP(o.p[k], kMass[k]);
                h0->Fill(o.bn0[k] - bp); h1->Fill(o.bn0[k] + o.db[k] - bp);
            }
            ax(h0, Form("#beta_{new} - #beta(p) for %s, %.0f < p < %.0f MeV/c", kLabel[k], plo[s], phi[s]), "particles (all rows)");
            h0->SetLineColor(kBlack); h0->SetLineWidth(2); h1->SetLineColor(kRed + 1); h1->SetLineWidth(2); h1->SetLineStyle(2);
            h0->SetMinimum(0.5); h0->SetMaximum(200 * std::max(h0->GetMaximum(), h1->GetMaximum()));
            h0->GetXaxis()->SetNdivisions(505);
            h0->Draw("hist"); h1->Draw("hist same");
            // polygon edges at the bin centre
            double lo = -1, hi = -1; for (double bb = 0.3; bb < 1.3; bb += 0.0001) { if (cut->IsInside(pm, bb)) { if (lo < 0) lo = bb; hi = bb; } }
            TLine e1(lo - b0, 0.5, lo - b0, h0->GetMaximum()); e1.SetLineColor(kGreen + 2); e1.SetLineWidth(2); e1.DrawClone();
            TLine e2(hi - b0, 0.5, hi - b0, h0->GetMaximum()); e2.SetLineColor(kGreen + 2); e2.SetLineWidth(2); e2.DrawClone();
            TLegend l(0.17, 0.70, 0.80, 0.88); l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextSize(0.04);
            l.AddEntry(h0, "stored (one segment)", "l"); l.AddEntry(h1, "two segments", "l"); l.AddEntry(&e1, Form("data polygon edges at p = %.0f", pm), "l");
            l.DrawClone();
            txt(0.14, 0.92, "DATA lep_034, PipPimEpEm_ID; vz > -160, PT3, isBest #geq 0", 0.042);
        }
        c.SaveAs(Form("%stofmom_slice.png", kFig));
        CF->Close();
    }
    printf("done\n");
    return 0;
}

}  // namespace tofmom

int sim_tofmom() { return tofmom::Main(); }
