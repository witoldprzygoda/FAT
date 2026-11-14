// ============================================================================
// sim_eta_enriched.C  --  data against simulation in eta-enriched selections below 1 GeV/c2
// ============================================================================
// User guidance (2026-11-12): SMASH describes the multi-pion background above
// the eta structures badly, so the simulation is NOT judged on the full
// spectrum. The comparison uses selections that keep mainly p p eta (protons
// not measured, all or almost all eta decay products measured) and some
// p p eta X, in M4 = M(pi+pi-e+e-) < 1 GeV/c2. A mismatch is a result for
// discussion, not a failure.
//
// Simulation: SMASH std (owned, mode A) + lep (owned, mode A) + eta (owned)
// e+e- rows, weight w (PDG BR), sigma = sum w / L_sim, L_sim = 221.05 nb^-1 per
// set; pure channel signals and physical background (templates 1 and 2), no
// combinatorial background. Ghost-flagged true pairs: variant (a)
// (sim_ghostvar.h). X17: shape for BR = 1e-5, information only, never summed.
// Data: (N - CB) / L_data, CB = 2 sqrt(N++ N--) (10 MeV bins in M4), L_data =
// 6080.5 nb^-1 UNCONFIRMED. All weighted sums in double.
//
// Selections (all: pre = isBest >= 0, PT3, vz > -160 mm, corrected momenta):
//   E1    slides:              oa_lab < 50, M(pipi) < 0.46, eta_rest > 140, OA > 2
//   E2    slides without OA:   oa_lab < 50, M(pipi) < 0.46, eta_rest > 140
//   E3    old code:            oa_lab < 50, M(pipi) < 0.42, eta_rest > 140
//   E3oa  old code + OA > 2    (OA is always shown with and without)
//   E4    eta band:            oa_lab < 50, eta_rest > 140
//   E4oa  eta band + OA > 2
//   E1mm, E2mm  E1, E2 + MM4 < cut (cut chosen here from the simulation, part 1)
//   E1v, E2v    E1, E2 + both pions in the vertex fit
//
// Physical background split: "eta-related" = a row with at least one track
// from an eta of the generator (geninfo 1717, inherited by every descendant;
// parent or grandparent GEANT id 17; phi -> gamma eta 550117; lep-set forced
// eta Dalitz 175101; eta-set forced decays), "multi-pion" = all other rows.
// The flags are computed from scan_nt (the compact files have no geninfo) and
// cached in /hdd3/.../scan/eta_enriched_etarel.root, aligned row by row with
// the compact trees (verified on pat_entry and w).
// Purity of "p p eta (X)" = (signal classes 2-5 + eta-related background) /
// (T1 + T2), event-group errors (simeng::Efficiency).
//
// Parts (bits): 1 MM4 study and choice of the MM4 cut (always computed; bit 1
// draws), 2 M4 comparison, regions, purity, 4 template fits, 8 subsystems.
// Output: figs/eta_enriched_*.png (pdf for the main ones), tables/eta_enriched_*.txt
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_enriched", kTRUE)' 'plots/eta_sim/sim_eta_enriched.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TTree.h"

#include "SimEngine.h"
#include "sim_ghostvar.h"

namespace seen {
using namespace simeng;

const char* const kCache = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/eta_enriched_etarel.root";
const double kLoadLo = 0.0, kLoadHi = 1.10;     // M4 range kept in memory (every selection uses M4 < 1.0)
const double kScaleK0S = 0.9952;                // sim M4 x 0.9952: K0S scale sim/data 1.0048 (review finding 8)
const double kBR17 = 1e-5;                      // X17 illustration
const char* const kLdTex = "L_{data} = 6080.5 nb^{-1} (UNCONFIRMED)";
const char* const kLdTxt = "L_data = 6080.5 nb^-1 (UNCONFIRMED)";
double gMMcut = -1;                             // MM4 cut, chosen in part 1

FILE* gTab = nullptr;                           // current table file (P writes to stdout and to it)
void P(const char* s) { fputs(s, stdout); if (gTab) fputs(s, gTab); }

// ---------------------------------------------------------------- selections
struct SelDef { const char* tag; const char* name; int base, oa, mm, vtx; };
const int kNSel = 10;
const SelDef kSel[kNSel] = {
    {"E1", "slides", 1, 1, 0, 0},
    {"E2", "slides without OA", 1, 0, 0, 0},
    {"E3", "old code", 2, 0, 0, 0},
    {"E3oa", "old code + OA > 2", 2, 1, 0, 0},
    {"E4", "eta band", 3, 0, 0, 0},
    {"E4oa", "eta band + OA > 2", 3, 1, 0, 0},
    {"E1mm", "E1 + MM4 cut", 1, 1, 1, 0},
    {"E2mm", "E2 + MM4 cut", 1, 0, 1, 0},
    {"E1v", "E1 + both pions in vertex fit", 1, 1, 0, 1},
    {"E2v", "E2 + both pions in vertex fit", 1, 0, 0, 1}};
enum { sE1 = 0, sE2, sE3, sE3oa, sE4, sE4oa, sE1mm, sE2mm, sE1v, sE2v };

template <class R> inline bool PassSel(int is, const R& r) {
    const SelDef& s = kSel[is];
    if (!PassPre(r)) return false;
    if (!(r.oa_lab < 50.f && r.oa_eta > 140.f)) return false;
    if (s.base == 1 && !(r.mpipi < 0.46f)) return false;
    if (s.base == 2 && !(r.mpipi < 0.42f)) return false;
    if (s.oa && !(r.oa > 2.f)) return false;
    if (s.mm && !(r.mm4 < (float)gMMcut)) return false;
    if (s.vtx && !(r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1)) return false;
    return true;
}
std::string SelLatex(int is) {
    const SelDef& s = kSel[is];
    std::string t = Form("%s (%s): #alpha_{lab} < 50#circ, #alpha_{#eta} > 140#circ", s.tag, s.name);
    if (s.base == 1) t += ", M_{#pi#pi} < 0.46";
    if (s.base == 2) t += ", M_{#pi#pi} < 0.42";
    if (s.oa) t += ", OA > 2#circ";
    if (s.mm) t += Form(", MM_{4} < %.2f", gMMcut);
    if (s.vtx) t += ", both #pi in vertex fit";
    return t;
}
std::string SelText(int is) {
    const SelDef& s = kSel[is];
    std::string t = "oa_lab < 50, eta_rest > 140";
    if (s.base == 1) t += ", M(pipi) < 0.46";
    if (s.base == 2) t += ", M(pipi) < 0.42";
    if (s.oa) t += ", OA > 2";
    if (s.mm) t += Form(", MM4 < %.2f", gMMcut);
    if (s.vtx) t += ", both pions in vertex fit";
    return t;
}
Sel SimSel(int is) { return [is](const Row& r) { return PassSel(is, r); }; }
scan::Sel DatSel(int is) { return [is](const scan::Row& r) { return PassSel(is, r); }; }
Sel SimSelM(int is, double lo, double hi) {
    return [is, lo, hi](const Row& r) { return r.m4 >= lo && r.m4 < hi && PassSel(is, r); };
}
scan::Sel DatSelM(int is, double lo, double hi) {
    return [is, lo, hi](const scan::Row& r) { return r.m4 >= lo && r.m4 < hi && PassSel(is, r); };
}

// ---------------------------------------------------------------- regions of M4
const int kNReg = 7;   // contiguous 0.35-1.00; "all" = their union
struct Reg { const char* tag; const char* name; double lo, hi; };
const Reg kReg[kNReg + 1] = {{"below", "below bump 0.35-0.455", 0.35, 0.455},
                             {"bump", "bump window 0.455-0.5225", scan::kBumpLo, scan::kBumpHi},
                             {"gap", "gap 0.5225-0.5336", scan::kBumpHi, scan::kPeakLo},
                             {"peak", "peak window 0.5336-0.5542", scan::kPeakLo, scan::kPeakHi},
                             {"post", "0.5542-0.56", scan::kPeakHi, 0.56},
                             {"r56", "0.56-0.70", 0.56, 0.70},
                             {"r70", "0.70-1.00", 0.70, 1.00},
                             {"all", "all 0.35-1.00", 0.35, 1.00}};
inline double RegX(double m) {
    if (!(m >= 0.35 && m < 1.00)) return -1.;
    for (int k = 0; k < kNReg; ++k) if (m >= kReg[k].lo && m < kReg[k].hi) return k + 0.5;
    return -1.;
}
Binning RegBinning(double scale = 1.) {
    Binning b; b.nbin = kNReg; b.lo = 0; b.hi = kNReg;
    b.x = [scale](const Row& r) { return RegX(r.m4 * scale); };
    return b;
}
// regions quoted in the report (index into kReg; 7 = all)
const int kNQ = 6;
const int kQ[kNQ] = {0, 1, 3, 5, 6, 7};

// ---------------------------------------------------------------- eta-related flags
// bits: 1 ep, 2 em, 4 pi+, 8 pi- track related to an eta; 16 a track of phi -> gamma eta (info)
std::vector<unsigned char> gER[3];
const SimData* gD = nullptr;
inline unsigned char ER(const Row& r) {
    if (r.set < 1 || r.set > 3) return 0;
    const std::vector<Row>& v = gD->t[TreeIndex(r.set, kUS)];
    const std::vector<unsigned char>& e = gER[r.set - 1];
    if (v.empty() || &r < v.data() || &r >= v.data() + v.size()) return 0;
    return e[(size_t)(&r - v.data())];
}
inline bool EtaRel(const Row& r) { return (ER(r) & 15) != 0; }
inline bool SigCls(const Row& r) { return r.tmpl == 1 && r.cls >= 2 && r.cls <= 5; }
inline bool EtaX(const Row& r) { return SigCls(r) || EtaRel(r); }   // "p p eta (X)"

inline bool EtaRelTrk(float g0, float par, float gp) {
    return g0 == 1717.f || g0 == 550117.f || g0 == 175101.f || g0 == simtruth::kEtaPPee || g0 == simtruth::kEtaX17 ||
           par == 17.f || gp == 17.f;
}

// Build the cache: one pass over scan_nt of the three e+e- outputs in parallel with the compact trees
bool BuildCache(FILE* log) {
    const char* nt[3] = {"output_std_pippimepem.root", "output_lep_pippimepem.root", "output_eta_pippimepem.root"};
    const char* cf[3] = {kCompact, kCompact, kCompactEta};
    TFile* fo = TFile::Open(kCache, "RECREATE");
    if (!fo || fo->IsZombie()) return false;
    bool ok = true;
    for (int s = 0; s < 3; ++s) {
        TStopwatch sw;
        TFile* fn = TFile::Open((std::string(kProdDir) + "/" + nt[s]).c_str(), "READ");
        TFile* fc = TFile::Open(cf[s], "READ");
        TTree* tn = fn ? (TTree*)fn->Get("scan_nt") : nullptr;
        TTree* tc = fc ? (TTree*)fc->Get(TreeName(TreeIndex(s + 1, kUS))) : nullptr;
        if (!tn || !tc) { ok = false; continue; }
        Float_t g[4][3], pe = 0, sw_ = 0;
        const char* sl[4] = {"ep", "em", "pip", "pim"};
        tn->SetBranchStatus("*", 0);
        for (int k = 0; k < 4; ++k) {
            const char* c[3] = {"_sim_geninfo", "_sim_parentid", "_sim_grandparentid"};
            for (int j = 0; j < 3; ++j) {
                const std::string b = std::string(sl[k]) + c[j];
                tn->SetBranchStatus(b.c_str(), 1);
                tn->SetBranchAddress(b.c_str(), &g[k][j]);
            }
        }
        tn->SetBranchStatus("pat_entry", 1); tn->SetBranchAddress("pat_entry", &pe);
        tn->SetBranchStatus("sim_w", 1); tn->SetBranchAddress("sim_w", &sw_);
        tn->SetCacheSize(200 * 1024 * 1024);
        Int_t cPat = 0, cA = 0, cB = 0; Float_t cM4 = 0; Double_t cW = 0;
        tc->SetBranchStatus("*", 0);
        for (const char* b : {"pat_entry", "ownA", "ownB", "m4", "w"}) tc->SetBranchStatus(b, 1);
        tc->SetBranchAddress("pat_entry", &cPat); tc->SetBranchAddress("ownA", &cA); tc->SetBranchAddress("ownB", &cB);
        tc->SetBranchAddress("m4", &cM4); tc->SetBranchAddress("w", &cW);
        const Long64_t n = tn->GetEntries();
        if (n != tc->GetEntries()) ok = false;
        fo->cd();
        TTree* to = new TTree(Form("er_%s", SetName(s + 1)), "eta-related bits per compact row (e+e-)");
        UChar_t er = 0; Float_t m4 = 0; Int_t own = 0, pat = 0;
        to->Branch("er", &er, "er/b"); to->Branch("m4", &m4, "m4/F");
        to->Branch("own", &own, "own/I"); to->Branch("pat", &pat, "pat/I");
        Long64_t bad = 0, nrel = 0;
        for (Long64_t i = 0; i < n; ++i) {
            tn->GetEntry(i); tc->GetEntry(i);
            if (cPat != (int)std::lround(pe) || cW != (double)sw_) ++bad;
            er = 0;
            for (int k = 0; k < 4; ++k) {
                if (EtaRelTrk(g[k][0], g[k][1], g[k][2])) er |= (1 << k);
                if (g[k][0] == 550117.f) er |= 16;
            }
            if (er & 15) ++nrel;
            m4 = cM4; own = (cA ? 1 : 0) | (cB ? 2 : 0); pat = cPat;
            to->Fill();
        }
        to->Write();
        const char* msg = Form("cache %s: %lld rows (compact %lld), mismatches pat_entry/w %lld, eta-related rows %lld, %.1f s\n",
                               SetName(s + 1), n, tc->GetEntries(), bad, nrel, sw.RealTime());
        printf("%s", msg); if (log) fputs(msg, log);
        if (bad) ok = false;
        fn->Close(); fc->Close(); delete fn; delete fc;
    }
    fo->Close(); delete fo;
    return ok;
}
// align the cached bits with the loaded rows (same filter as simeng::Load)
bool AttachER(const SimData& d, FILE* log) {
    TFile* f = TFile::Open(kCache, "READ");
    if (!f || f->IsZombie()) return false;
    bool ok = true;
    for (int s = 0; s < 3; ++s) {
        TTree* t = (TTree*)f->Get(Form("er_%s", SetName(s + 1)));
        if (!t) { ok = false; continue; }
        UChar_t er = 0; Float_t m4 = 0; Int_t own = 0, pat = 0;
        t->SetBranchAddress("er", &er); t->SetBranchAddress("m4", &m4);
        t->SetBranchAddress("own", &own); t->SetBranchAddress("pat", &pat);
        const std::vector<Row>& v = d.t[TreeIndex(s + 1, kUS)];
        gER[s].clear(); gER[s].reserve(v.size());
        Long64_t bad = 0;
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (!(m4 >= kLoadLo && m4 < kLoadHi)) continue;
            if (!own) continue;
            const size_t j = gER[s].size();
            if (j >= v.size() || v[j].pat_entry != pat) ++bad;
            gER[s].push_back(er);
        }
        const char* msg = Form("attach %s: %zu cached rows for %zu loaded rows, pat_entry mismatches %lld\n", SetName(s + 1),
                               gER[s].size(), v.size(), bad);
        printf("%s", msg); if (log) fputs(msg, log);
        if (bad || gER[s].size() != v.size()) ok = false;
    }
    f->Close(); delete f;
    return ok;
}

// ---------------------------------------------------------------- components
enum { cBD = 0, cBC, cPC, cPD, cKE, cKM, cSum, cGE, cGM, cSumA, cX17, cF, cDf, cKEpi, cKEee, cKEoth, cN };
std::vector<Comp> Comps() {
    std::vector<Comp> c(cN);
    c[cBD] = CompBumpDalitz(); c[cBC] = CompBumpConv(); c[cPC] = CompPeakConv(); c[cPD] = CompPeakDirect();
    c[cKE] = Comp("phys. bkg, eta-related", [](const Row& r) { return r.tmpl == 2 && EtaRel(r); });
    c[cKM] = Comp("phys. bkg, multi-pion", [](const Row& r) { return r.tmpl == 2 && !EtaRel(r); });
    c[cSum] = CompSim();
    c[cGE] = Comp("ghost-flagged true pairs, eta-related", [](const Row& r) { return GhostPair(r) && EtaRel(r); });
    c[cGM] = Comp("ghost-flagged true pairs, multi-pion", [](const Row& r) { return GhostPair(r) && !EtaRel(r); });
    c[cSumA] = CompSimVar(1);
    c[cX17] = CompX17(kBR17);
    c[cF] = Comp("info: fake lepton (T3 without ghost pairs)", [](const Row& r) { return r.tmpl == 3 && !GhostPair(r); });
    c[cDf] = Comp("info: different sources (T4)", [](const Row& r) { return r.tmpl == 4; });
    c[cKEpi] = Comp("  eta-rel.: both pions from an eta", [](const Row& r) { return r.tmpl == 2 && (ER(r) & 12) == 12; });
    c[cKEee] = Comp("  eta-rel.: both leptons from an eta (not both pions)",
                    [](const Row& r) { return r.tmpl == 2 && (ER(r) & 12) != 12 && (ER(r) & 3) == 3; });
    c[cKEoth] = Comp("  eta-rel.: other (one or two eta tracks)",
                     [](const Row& r) { return r.tmpl == 2 && (ER(r) & 12) != 12 && (ER(r) & 3) != 3 && (ER(r) & 15) != 0; });
    return c;
}
const char* CompShort(int c) {
    static const char* n[cN] = {"bump Dalitz", "bump conv.", "peak conv.", "peak direct", "bkg eta-rel.", "bkg multi-pion",
                                "sim sum T1+T2", "ghost eta-rel.", "ghost multi-pion", "sum (a)", "X17 1e-5", "fake lepton",
                                "diff. sources", "eta-rel: pions", "eta-rel: leptons", "eta-rel: other"};
    return n[c];
}
// fit components (combined so that their event-group variance is exact)
enum { fBump = 0, fPeak, fPD, fPC, fK, fKE, fKM, fKA, fKEA, fKMA, fN };
std::vector<Comp> FitComps() {
    std::vector<Comp> c(fN);
    c[fBump] = Comp("bump (Dalitz + conv.)", [](const Row& r) { return r.tmpl == 1 && (r.cls == 2 || r.cls == 3); });
    c[fPeak] = CompPeakTotal();
    c[fPD] = CompPeakDirect();
    c[fPC] = CompPeakConv();
    c[fK] = CompPhysBkg();
    c[fKE] = Comp("bkg eta-rel.", [](const Row& r) { return r.tmpl == 2 && EtaRel(r); });
    c[fKM] = Comp("bkg multi-pion", [](const Row& r) { return r.tmpl == 2 && !EtaRel(r); });
    c[fKA] = Comp("bkg (a)", [](const Row& r) { return r.tmpl == 2 || GhostPair(r); });
    c[fKEA] = Comp("bkg eta-rel. (a)", [](const Row& r) { return (r.tmpl == 2 || GhostPair(r)) && EtaRel(r); });
    c[fKMA] = Comp("bkg multi-pion (a)", [](const Row& r) { return (r.tmpl == 2 || GhostPair(r)) && !EtaRel(r); });
    return c;
}

// ---------------------------------------------------------------- data histogram of any variable
// N per bin with the event-group variance, N++ and N-- per bin, CB = 2 sqrt(N++ N--) per bin
struct DH {
    int nb = 0; double lo = 0, hi = 0, bw = 0;
    std::vector<double> n, vn, pp, mm;
    double cb(int i) const { return 2. * std::sqrt(pp[i] * mm[i]); }
    double sig(int i) const { return n[i] - cb(i); }
    double sigVar(int i) const { return vn[i] + pp[i] + mm[i]; }
    double Sum(bool sg) const { double s = 0; for (int i = 0; i < nb; ++i) s += sg ? sig(i) : n[i]; return s; }
    double SumVar() const { double s = 0; for (int i = 0; i < nb; ++i) s += sigVar(i); return s; }
    // unit 2: rows / L_data / bw ; what 0 N, 1 CB, 2 N - CB
    TH1D* Hist(const char* name, int what, double scale) const {
        TH1D* h = new TH1D(name, "", nb, lo, hi); h->SetDirectory(nullptr);
        for (int i = 0; i < nb; ++i) {
            const double v = what == 0 ? n[i] : what == 1 ? cb(i) : sig(i);
            const double e = what == 0 ? std::sqrt(vn[i]) : what == 1 ? std::sqrt(pp[i] + mm[i]) : std::sqrt(sigVar(i));
            h->SetBinContent(i + 1, v * scale); h->SetBinError(i + 1, e * scale);
        }
        return h;
    }
};
typedef std::function<double(const scan::Row&)> DVar;
// several variables of one selection in one pass
std::vector<DH> FillDH(const scan::Data& dd, const scan::Sel& sel, const std::vector<DVar>& xs,
                       const std::vector<Binning>& bs) {
    const int nv = (int)xs.size();
    std::vector<DH> o(nv);
    for (int k = 0; k < nv; ++k) {
        o[k].nb = bs[k].nbin; o[k].lo = bs[k].lo; o[k].hi = bs[k].hi; o[k].bw = bs[k].bw();
        o[k].n.assign(o[k].nb, 0.); o[k].vn.assign(o[k].nb, 0.); o[k].pp.assign(o[k].nb, 0.); o[k].mm.assign(o[k].nb, 0.);
    }
    for (int smp = 0; smp < 3; ++smp) {
        const std::vector<scan::Row>& v = dd.s[smp];
        std::vector<std::vector<std::pair<int, double>>> g(nv);
        int cur = -1;
        auto flush = [&]() {
            for (int k = 0; k < nv; ++k) {
                if (g[k].empty()) continue;
                std::sort(g[k].begin(), g[k].end());
                size_t i = 0;
                while (i < g[k].size()) {
                    const int b = g[k][i].first; double a = 0;
                    while (i < g[k].size() && g[k][i].first == b) { a += g[k][i].second; ++i; }
                    o[k].vn[b] += a * a;
                }
                g[k].clear();
            }
        };
        for (const scan::Row& r : v) {
            if (smp == 0 && r.grp != cur) { flush(); cur = r.grp; }
            if (!r.cor_ok) continue;
            if (!sel(r)) continue;
            for (int k = 0; k < nv; ++k) {
                const double x = xs[k](r);
                if (!(x >= o[k].lo && x < o[k].hi)) continue;
                const int b = std::min(o[k].nb - 1, (int)((x - o[k].lo) / o[k].bw));
                if (smp == 0) { o[k].n[b] += 1; g[k].push_back(std::make_pair(b, 1.)); }
                else if (smp == 1) o[k].pp[b] += 1;
                else o[k].mm[b] += 1;
            }
        }
        if (smp == 0) flush();
    }
    return o;
}

// ---------------------------------------------------------------- plotting helpers
void StyleH(TH1D* h, int col, int lw, int ls, int fill = 0) {
    h->SetLineColor(col); h->SetLineWidth(lw); h->SetLineStyle(ls);
    h->SetMarkerColor(col); h->SetMarkerSize(0);
    if (fill) { h->SetFillColor(fill); h->SetFillStyle(1001); } else h->SetFillStyle(0);
}
void StyleM(TH1D* h, int col, int mst, double ms = 1.0) {
    h->SetMarkerStyle(mst); h->SetMarkerColor(col); h->SetLineColor(col); h->SetMarkerSize(ms); h->SetLineWidth(2);
}
void Axes(TH1D* h, const char* xt, const char* yt, double ts = 0.055, double yoff = 1.3) {
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(ts); h->GetYaxis()->SetTitleSize(ts);
    h->GetXaxis()->SetLabelSize(ts * 0.87); h->GetYaxis()->SetLabelSize(ts * 0.87);
    h->GetXaxis()->SetTitleOffset(1.0); h->GetYaxis()->SetTitleOffset(yoff);
    h->GetYaxis()->SetMaxDigits(3);
    h->GetXaxis()->SetNdivisions(506);
}
void Txt(double x, double y, const char* s, double size, int col = 1, int align = 11) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(align); t.DrawLatex(x, y, s);
}
TPad* MkPad(const char* name, double x1, double y1, double x2, double y2, bool logy, double l = 0.14, double r = 0.03,
            double t = 0.08, double b = 0.03) {
    TPad* p = new TPad(name, "", x1, y1, x2, y2);
    p->SetLeftMargin(l); p->SetRightMargin(r); p->SetTopMargin(t); p->SetBottomMargin(b);
    p->SetLogy(logy);
    p->Draw(); p->cd();
    return p;
}
void VLines(const std::vector<double>& xs, double ylo, double yhi, int col = kGray + 1, int st = 3) {
    for (double x : xs) { TLine* l = new TLine(x, ylo, x, yhi); l->SetLineStyle(st); l->SetLineColor(col); l->Draw(); }
}
TH1D* RatioH(const TH1D* a, const TH1D* b, const char* name) {
    TH1D* r = (TH1D*)a->Clone(name); r->Reset(); r->SetDirectory(nullptr);
    for (int i = 1; i <= a->GetNbinsX(); ++i) {
        const double x = a->GetBinContent(i), y = b->GetBinContent(i), ex = a->GetBinError(i), ey = b->GetBinError(i);
        if (!(y > 0)) continue;
        const double v = x / y;
        r->SetBinContent(i, v);
        r->SetBinError(i, std::sqrt(ex * ex / (y * y) + v * v * ey * ey / (y * y)));
    }
    return r;
}
double HMax(const TH1D* h, double lo = -1e30, double hi = 1e30, bool err = true) {
    double m = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        const double x = h->GetBinCenter(i);
        if (x < lo || x > hi) continue;
        m = std::max(m, h->GetBinContent(i) + (err ? h->GetBinError(i) : 0.));
    }
    return m;
}
double HMinPos(const TH1D* h) {
    double m = 1e30;
    for (int i = 1; i <= h->GetNbinsX(); ++i) if (h->GetBinContent(i) > 0) m = std::min(m, h->GetBinContent(i));
    return m;
}
TH1D* Dummy(const char* n, int lc, int lw, int ls, int fc, int mst = 0, int mc = 1) {
    TH1D* h = new TH1D(n, "", 1, 0, 1); h->SetDirectory(nullptr);
    h->SetLineColor(lc); h->SetLineWidth(lw); h->SetLineStyle(ls);
    if (fc) { h->SetFillColor(fc); h->SetFillStyle(1001); } else h->SetFillStyle(0);
    if (mst) { h->SetMarkerStyle(mst); h->SetMarkerColor(mc); h->SetMarkerSize(1.6); }
    return h;
}
const int kColData = kBlack, kColCB = kRed + 1, kColSig = kBlue + 1;
// stacking order of the simulated components (bottom to top) and their fill colours
const int kStack[6] = {cKM, cKE, cBD, cBC, cPC, cPD};
const int kStackLine[6] = {kColPhysBkg, kColEtaBkg, kColBumpDalitz, kColBumpConv, kColPeakConv, kColPeakDirect};
const int kStackFill[6] = {kColPhysFill, kColEtaFill, kColBumpDalitz, kColBumpConv, kColPeakConv, kColPeakDirect};
const char* kStackLeg[6] = {"phys. bkg, multi-pion", "phys. bkg, eta-related", "bump: #pi^{0} Dalitz",
                            "bump: conversion", "peak: conversion (#pi^{+}#pi^{-}#gamma)", "peak: direct #pi^{+}#pi^{-}e^{+}e^{-} (M_{ee} not QED)"};
void SharedLegend(double x1, double y1, double x2, double y2, int ncol, double ts, bool withData = true, bool withCB = true,
                  bool withX17 = true) {
    TLegend* lg = new TLegend(x1, y1, x2, y2);
    lg->SetNColumns(ncol); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(ts);
    if (withData) {
        if (withCB) lg->AddEntry(Dummy("lgN", kColData, 2, 1, 0, 20, kColData), "data N", "p");
        if (withCB) lg->AddEntry(Dummy("lgCB", kColCB, 2, 1, 0, 21, kColCB), "data CB = 2#sqrt{N_{++}N_{--}}", "p");
        lg->AddEntry(Dummy("lgS", kColSig, 2, 1, 0, 20, kColSig), "data N - CB", "p");
    }
    for (int k = 5; k >= 0; --k) lg->AddEntry(Dummy(Form("lgK%d", k), kStackLine[k], 2, 1, kStackFill[k]), kStackLeg[k], "f");
    lg->AddEntry(Dummy("lgSum", kColSim, 4, 1, 0), "simulation T1 + T2 (sum)", "l");
    lg->AddEntry(Dummy("lgSumA", kColSim, 3, 2, 0), "variant (a): + ghost-flagged true pairs", "l");
    if (withX17) lg->AddEntry(Dummy("lgX", kColX17, 3, 2, 0), "X17 shape, BR 10^{-5} (information, not in the sum)", "l");
    lg->Draw();
}
const char* kXM4 = "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]";
const char* kYM4 = "d#sigma/dM [nb/(GeV/c^{2})]";
std::string NormLine() {
    return Form("data: rows / L_{data} / bin width, %s; simulation: SMASH std + lep + eta, owned rows (mode A), w (PDG), "
                "#Sigma w / L_{sim} / bin width, L_{sim} = %.2f nb^{-1} per set; no efficiency correction; no CB in simulation",
                kLdTex, kLsim);
}
const char* PreLine() { return "all selections: pre = isBest #geq 0, PT3, v_{z} > -160 mm, corrected momenta (cor_ok), e^{+}e^{-}"; }

// ---------------------------------------------------------------- template fit
struct FitRes {
    bool ok = false; int np = 0, npts = 0, ndf = 0, iter = 0;
    double chi2 = 0, a[6], ea[6], cov[6][6];
};
// y: data N - CB per bin (rows); sp: data spectra for the covariance (Poisson N + interpolated CB)
// T[k], TV[k]: template k in expected data rows and its variance per bin; bins: indices used
FitRes TemplateFit(const scan::Spectra& sp, const std::vector<int>& bins, const std::vector<std::vector<double>>& T,
                   const std::vector<std::vector<double>>& TV) {
    FitRes o;
    const int np = (int)T.size(), n = (int)bins.size();
    o.np = np; o.npts = n; o.ndf = n - np;
    for (int k = 0; k < 6; ++k) { o.a[k] = 1; o.ea[k] = 0; for (int l = 0; l < 6; ++l) o.cov[k][l] = 0; }
    std::vector<double> y(n), F(n * np), V0(n * n, 0.);
    for (int i = 0; i < n; ++i) {
        const int b = bins[i];
        y[i] = sp.sig[b];
        for (int k = 0; k < np; ++k) F[i * np + k] = T[k][b];
        for (int j = 0; j < n; ++j) V0[i * n + j] = sp.CbCov(b, bins[j]);
        V0[i * n + i] += std::max(1., sp.n[0][b]);
    }
    double a[6] = {1, 1, 1, 1, 1, 1};
    scan::LinRes r;
    for (int it = 0; it < 50; ++it) {
        std::vector<double> V = V0;
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < np; ++k) V[i * n + i] += a[k] * a[k] * TV[k][bins[i]];
        r = scan::LinFitCov(n, np, y.data(), V, F.data());
        if (!r.ok) return o;
        double d = 0;
        for (int k = 0; k < np; ++k) { d = std::max(d, std::fabs(r.par[k] - a[k])); a[k] = r.par[k]; }
        o.iter = it + 1;
        if (d < 1e-7) break;
    }
    std::vector<double> V = V0;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < np; ++k) V[i * n + i] += a[k] * a[k] * TV[k][bins[i]];
    r = scan::LinFitCov(n, np, y.data(), V, F.data());
    if (!r.ok) return o;
    o.ok = true; o.chi2 = r.chi2;
    for (int k = 0; k < np; ++k) {
        o.a[k] = r.par[k]; o.ea[k] = std::sqrt(std::max(0., r.cov[k][k]));
        for (int l = 0; l < np; ++l) o.cov[k][l] = r.cov[k][l];
    }
    return o;
}

// P_main by hand (the data-study peak fit): Gaussian with frozen mean and width + pol2 in the frozen range
// 0.5225-0.61 (bin centres), weighted linear least squares with the given variances. Used on the simulation
// (variances = template variances) and, as a cross-check, on the data (variances = N + CB, diagonal).
struct PM { bool ok = false; double S = 0, eS = 0, chi2 = 0; int ndf = 0; };
PM PmainLin(const std::vector<double>& y, const std::vector<double>& var, double lo, double bw, double mean, double sigma) {
    const double flo = scan::kFitLo, fhi = scan::kFitHi, x0 = 0.5 * (flo + fhi), dx = 0.5 * (fhi - flo);
    double vmax = 0;
    for (double v : var) vmax = std::max(vmax, v);
    std::vector<double> yy, vv, F;
    for (size_t i = 0; i < y.size(); ++i) {
        const double xc = lo + (i + 0.5) * bw;
        if (!(xc >= flo && xc < fhi)) continue;
        yy.push_back(y[i]); vv.push_back(std::max(var[i], 1e-6 * vmax + 1e-12));
        const double t = (xc - x0) / dx;
        F.push_back(scan::GaussFrac(xc - 0.5 * bw, xc + 0.5 * bw, mean, sigma));
        F.push_back(1.); F.push_back(t); F.push_back(t * t);
    }
    const scan::LinRes r = scan::LinFit((int)yy.size(), 4, yy.data(), vv.data(), F.data());
    PM o; o.ok = r.ok; o.S = r.par[0]; o.eS = std::sqrt(std::max(0., r.cov[0][0])); o.chi2 = r.chi2; o.ndf = (int)yy.size() - 4;
    return o;
}

}  // namespace seen

// ============================================================================
void sim_eta_enriched(int parts = 15) {
    using namespace simeng;
    using namespace seen;
    TStopwatch swAll;
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(0);
    TH1::AddDirectory(kFALSE);
    const std::string tab = kTabDir, fig = kFigDir;

    // ---------------- checks table ----------------
    FILE* fchk = fopen((tab + "/eta_enriched_checks.txt").c_str(), "w");
    fprintf(fchk, "# eta_enriched_checks: inputs, cache alignment and consistency checks of sim_eta_enriched.C\n");
    {
        TFile* f = TFile::Open(kCache, "READ");
        const bool have = f && !f->IsZombie() && f->Get("er_std") && f->Get("er_lep") && f->Get("er_eta");
        if (f) { f->Close(); delete f; }
        if (!have) {
            printf("building the eta-related cache %s\n", kCache);
            if (!BuildCache(fchk)) { printf("cache build FAILED\n"); fclose(fchk); return; }
        } else fprintf(fchk, "cache %s present (built by an earlier run of this macro)\n", kCache);
    }

    // ---------------- load ----------------
    SimData d;
    if (!Load(d, kCompact, kMaskUS, kLoadLo, kLoadHi, 1)) return;
    gD = &d;
    if (!AttachER(d, fchk)) { printf("eta-related flags not aligned: STOP\n"); fclose(fchk); return; }
    scan::Data dd;
    if (!scan::Load(dd, scan::kCompact, kLoadLo, kLoadHi)) return;
    SimOpt opt;   // std + lep + eta, mode A, w, e+e-
    const std::vector<Comp> comps = Comps();

    // consistency: E1 = preset slides, E3 = preset old code (data and simulation rows)
    {
        long long dif[2] = {0, 0}, nsel[2] = {0, 0};
        for (int k = 0; k < 3; ++k)
            for (const scan::Row& r : dd.s[k]) {
                if (PassSel(sE1, r) != PassPreset(kSlides, r)) dif[0]++;
                if (PassSel(sE3, r) != PassPreset(kOldCode, r)) dif[1]++;
                if (PassSel(sE1, r)) nsel[0]++;
                if (PassSel(sE3, r)) nsel[1]++;
            }
        long long sdif = 0;
        for (int k : {kStdUS, kLepUS, kEtaUS})
            for (const Row& r : d.t[k]) {
                if (PassSel(sE1, r) != PassPreset(kSlides, r)) sdif++;
                if (PassSel(sE3, r) != PassPreset(kOldCode, r)) sdif++;
            }
        fprintf(fchk, "E1 vs preset slides, E3 vs preset old code: data rows differing %lld / %lld (selected %lld / %lld); sim rows differing %lld\n",
                dif[0], dif[1], nsel[0], nsel[1], sdif);
        // eta-related share per set (owned rows, pre)
        for (int s = 1; s <= 3; ++s) {
            double wa = 0, we = 0;
            for (const Row& r : d.tree(s, kUS)) {
                if (!Owned(r, opt) || !PassPre(r) || !(r.tmpl == 2)) continue;
                wa += r.w; if (EtaRel(r)) we += r.w;
            }
            fprintf(fchk, "set %s: physical background (T2, pre, M4 < 1.1): sum w %.6g, eta-related share %.4f\n", SetName(s), wa,
                    wa > 0 ? we / wa : 0.);
        }
    }

    // ======================================================================
    // part 1: MM4 distributions and the choice of the MM4 cut
    // ======================================================================
    {
        TStopwatch sw;
        // fine scan binning (10 MeV) over the whole MM4 range; display binning 50 MeV in 1.60-3.30
        Binning bScan; bScan.nbin = 700; bScan.lo = -3.5; bScan.hi = 3.5; bScan.x = [](const Row& r) { return (double)r.mm4; };
        Binning bDisp; bDisp.nbin = 34; bDisp.lo = 1.60; bDisp.hi = 3.30; bDisp.x = bScan.x;
        std::vector<Comp> mc = {Comp("p p eta (X): signal classes 2-5 + eta-related bkg", [](const Row& r) { return EtaX(r) && (r.tmpl == 1 || r.tmpl == 2) && !X17Row(r); }),
                                Comp("signal classes 2-5", [](const Row& r) { return SigCls(r); }),
                                comps[cKE], comps[cKM], comps[cSum], comps[cBD], comps[cBC], comps[cPC], comps[cPD]};
        enum { mEX = 0, mSig, mKE, mKM, mSum, mBD, mBC, mPC, mPD, mN };
        const int ms[2] = {sE2, sE1};
        std::vector<WSpec> scanS[2], dispS[2];
        std::vector<DH> dScan[2], dDisp[2];
        for (int j = 0; j < 2; ++j) {
            const Sel sl = SimSelM(ms[j], 0.35, 1.00);
            scanS[j] = FillMany(d, sl, mc, opt, bScan);
            dispS[j] = FillMany(d, sl, mc, opt, bDisp);
            DVar xm = [](const scan::Row& r) { return (double)r.mm4; };
            dScan[j] = FillDH(dd, DatSelM(ms[j], 0.35, 1.00), {xm}, {bScan});
            dDisp[j] = FillDH(dd, DatSelM(ms[j], 0.35, 1.00), {xm}, {bDisp});
        }
        // cumulative fractions for cuts on a 0.05 grid
        auto cumS = [&](const WSpec& s, double c) {
            double a = 0; for (int i = 0; i < s.nbin; ++i) if (s.lo + (i + 1) * s.bw <= c + 1e-9) a += s.s[i];
            return s.S > 0 ? a / s.S : 0.;
        };
        auto cumD = [&](const DH& h, double c) {
            double a = 0, t = 0;
            for (int i = 0; i < h.nb; ++i) { t += h.sig(i); if (h.lo + (i + 1) * h.bw <= c + 1e-9) a += h.sig(i); }
            return t > 0 ? a / t : 0.;
        };
        gTab = fopen((tab + "/eta_enriched_mm4.txt").c_str(), "w");
        P("# eta_enriched_mm4: MM4 = missing mass of pi+pi-e+e- [GeV/c2], selections E2 and E1, M4 in 0.35-1.00\n");
        P("# fraction of each component with MM4 < c (simulation: sum w, std + lep + eta owned rows, mode A, w; data: N - CB with CB per 10 MeV MM4 bin)\n");
        P("# p p eta (X) = signal classes 2-5 + physical background rows with a track from an eta; multi-pion = other T2 rows\n");
        P("# exclusive p p eta: MM4 = M(pp) in 1.876-2.919 (sqrt(s) - m_eta); any row: MM4 <= sqrt(s) - M4\n");
        P("# RULE (decided here): the MM4 cut is the smallest value on a 0.05 grid that keeps >= 90 % of the simulated signal classes 2-5\n");
        P("#   (eta -> pi+pi-pi0 and eta -> pi+pi-gamma / pi+pi-e+e-, all products but a photon measured) in selection E2, M4 0.35-1.00\n");
        double cutRule = -1;
        for (int j = 0; j < 2; ++j) {
            P(Form("\n## selection %s (%s), M4 0.35-1.00\n", kSel[ms[j]].tag, SelText(ms[j]).c_str()));
            P(Form("%6s %10s %10s %10s %10s %10s %10s %12s\n", "cut c", "ppeta(X)", "sig 2-5", "bkg eta", "multi-pi", "sim sum", "data N-CB",
                   "purity(sim)"));
            for (double c = 2.30; c < 3.251; c += 0.05) {
                const double eEX = cumS(scanS[j][mEX], c), eSig = cumS(scanS[j][mSig], c), eKE = cumS(scanS[j][mKE], c),
                             eKM = cumS(scanS[j][mKM], c), eSum = cumS(scanS[j][mSum], c), eD = cumD(dScan[j][0], c);
                const double pur = (eSum * scanS[j][mSum].S) > 0 ? eEX * scanS[j][mEX].S / (eSum * scanS[j][mSum].S) : 0;
                P(Form("%6.2f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %12.4f\n", c, eEX, eSig, eKE, eKM, eSum, eD, pur));
                if (j == 0 && cutRule < 0 && eSig >= 0.90) cutRule = c;
            }
            P(Form("no cut: purity (sim) %.4f; sum w: ppeta(X) %.6g, multi-pion %.6g\n", scanS[j][mEX].S / scanS[j][mSum].S,
                   scanS[j][mEX].S, scanS[j][mKM].S));
        }
        if (cutRule < 0) cutRule = 3.25;   // never reached: effectively no cut
        gMMcut = std::floor(cutRule * 100 + 0.5) / 100.;
        P(Form("\nCHOSEN MM4 cut: MM4 < %.2f GeV/c2 (rule above, from the simulation only)\n", gMMcut));
        fclose(gTab); gTab = nullptr;

        if (parts & 1) {
            TCanvas* cv = new TCanvas("cMM", "", 2700, 1900);
            cv->cd();
            TPad* ph = MkPad("pHmm", 0, 0.88, 1, 1, false, 0, 0, 0, 0);
            Txt(0.01, 0.80, Form("Missing mass MM_{4} of #pi^{+}#pi^{-}e^{+}e^{-} in the eta-enriched selections E2 and E1, M_{4} = 0.35-1.00 GeV/c^{2}; chosen cut MM_{4} < %.2f GeV/c^{2}", gMMcut), 0.15);
            Txt(0.01, 0.60, NormLine().c_str(), 0.095);
            Txt(0.01, 0.43, Form("%s; exclusive p p #eta: MM_{4} = M_{pp} #in [1.876, 2.919]; rule: smallest cut keeping #geq 90 %% of the simulated signal classes 2-5 in E2", PreLine()), 0.095);
            TLegend* lg = new TLegend(0.0, 0.0, 1.0, 0.33);
            lg->SetNColumns(5); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.09);
            lg->AddEntry(Dummy("mN", kColData, 2, 1, 0, 20, kColData), "data N", "p");
            lg->AddEntry(Dummy("mCB", kColCB, 2, 1, 0, 21, kColCB), "data CB", "p");
            lg->AddEntry(Dummy("mS", kColSig, 2, 1, 0, 20, kColSig), "data N - CB", "p");
            lg->AddEntry(Dummy("mSum", kColSim, 4, 1, 0), "sim T1 + T2", "l");
            lg->AddEntry(Dummy("mEX", kGreen + 2, 4, 1, 0), "sim p p #eta (X): signal 2-5 + eta-related bkg", "l");
            lg->AddEntry(Dummy("mSig", kColBumpDalitz, 3, 1, 0), "sim signal classes 2-5", "l");
            lg->AddEntry(Dummy("mKE", kColEtaBkg, 3, 1, 0), "sim bkg eta-related", "l");
            lg->AddEntry(Dummy("mKM", kColPhysBkg, 3, 1, kColPhysFill), "sim bkg multi-pion", "lf");
            lg->AddEntry(Dummy("mPur", kMagenta + 2, 3, 7, 0), "sim purity p p #eta (X) for MM_{4} < c", "l");
            lg->AddEntry(Dummy("mCut", kRed, 3, 2, 0), "chosen cut", "l");
            lg->Draw();
            const char* yAbs = "d#sigma/dMM_{4} [nb/(GeV/c^{2})]";
            // common y range of the absolute panels
            double yhi = 0, ylo = 1e30;
            std::vector<TH1D*> hd[2], hs[2];
            for (int j = 0; j < 2; ++j) {
                const double sc = 1. / (kLdata * bDisp.bw());
                hd[j] = {dDisp[j][0].Hist(Form("mmN%d", j), 0, sc), dDisp[j][0].Hist(Form("mmC%d", j), 1, sc),
                         dDisp[j][0].Hist(Form("mmS%d", j), 2, sc)};
                for (int c = 0; c < mN; ++c) hs[j].push_back(dispS[j][c].Hist(Form("mmH%d_%d", j, c), 2));
                yhi = std::max(yhi, HMax(hd[j][0]));
                ylo = std::min(ylo, HMinPos(hs[j][mSig]));
            }
            for (int j = 0; j < 2; ++j) {
                const double y1 = j == 0 ? 0.44 : 0.0, y2 = j == 0 ? 0.88 : 0.44;
                // (a) absolute, log
                cv->cd(); MkPad(Form("pMa%d", j), 0.0, y1, 0.34, y2, true, 0.14, 0.03, 0.10, 0.13);
                TH1D* fr = (TH1D*)hs[j][mSum]->Clone(Form("mfr%d", j)); fr->Reset(); Axes(fr, "MM_{4} [GeV/c^{2}]", yAbs, 0.05);
                fr->SetMinimum(std::max(ylo * 0.3, yhi * 1e-5)); fr->SetMaximum(yhi * 5); fr->Draw("axis");
                StyleH(hs[j][mKM], kColPhysBkg, 2, 1, kColPhysFill); hs[j][mKM]->Draw("hist same");
                StyleH(hs[j][mKE], kColEtaBkg, 3, 1); hs[j][mKE]->Draw("hist same");
                StyleH(hs[j][mSig], kColBumpDalitz, 3, 1); hs[j][mSig]->Draw("hist same");
                StyleH(hs[j][mEX], kGreen + 2, 4, 1); hs[j][mEX]->Draw("hist same");
                StyleH(hs[j][mSum], kColSim, 4, 1); hs[j][mSum]->Draw("hist same");
                StyleM(hd[j][0], kColData, 20); hd[j][0]->Draw("E1 same");
                StyleM(hd[j][1], kColCB, 21); hd[j][1]->Draw("E1 same");
                StyleM(hd[j][2], kColSig, 20); hd[j][2]->Draw("E1 same");
                VLines({1.876, 2.919}, fr->GetMinimum(), yhi * 5, kGray + 2, 3);
                VLines({gMMcut}, fr->GetMinimum(), yhi * 5, kRed, 2);
                fr->Draw("axis same");
                Txt(0.15, 0.92, Form("(%s) %s: absolute, log; 50 MeV bins", j == 0 ? "a" : "d", kSel[ms[j]].tag), 0.05);
                Txt(0.15, 0.86, SelLatex(ms[j]).c_str(), 0.030);
                Txt(0.17, 0.80, "dotted: 1.876 and 2.919 (exclusive p p #eta)", 0.033, kGray + 2);
                // (b) normalised shapes
                cv->cd(); MkPad(Form("pMb%d", j), 0.34, y1, 0.67, y2, false, 0.14, 0.03, 0.10, 0.13);
                std::vector<TH1D*> hn;
                for (TH1D* h : {hd[j][2], hs[j][mEX], hs[j][mKM], hs[j][mSum]}) {
                    TH1D* c = (TH1D*)h->Clone(Form("%s_n", h->GetName()));
                    const double in = c->Integral();
                    if (in > 0) c->Scale(1. / in);
                    hn.push_back(c);
                }
                double ym = 0; for (TH1D* h : hn) ym = std::max(ym, HMax(h));
                TH1D* fb = (TH1D*)hn[0]->Clone(Form("mfb%d", j)); fb->Reset();
                Axes(fb, "MM_{4} [GeV/c^{2}]", "fraction per 50 MeV", 0.05);
                fb->SetMinimum(0); fb->SetMaximum(ym * 1.25); fb->Draw("axis");
                StyleH(hn[2], kColPhysBkg, 3, 1, 0); hn[2]->Draw("hist same");
                StyleH(hn[1], kGreen + 2, 4, 1); hn[1]->Draw("hist same");
                StyleH(hn[3], kColSim, 3, 2); hn[3]->Draw("hist same");
                StyleM(hn[0], kColSig, 20); hn[0]->Draw("E1 same");
                VLines({gMMcut}, 0, ym * 1.25, kRed, 2);
                fb->Draw("axis same");
                Txt(0.15, 0.92, Form("(%s) %s: shapes, each normalised to 1", j == 0 ? "b" : "e", kSel[ms[j]].tag), 0.05);
                // (c) cumulative fractions
                cv->cd(); MkPad(Form("pMc%d", j), 0.67, y1, 1.0, y2, false, 0.14, 0.03, 0.10, 0.13);
                TH1D* fc = new TH1D(Form("mfc%d", j), "", 10, 2.30, 3.25); fc->SetDirectory(nullptr);
                Axes(fc, "cut c: MM_{4} < c [GeV/c^{2}]", "fraction kept (M_{4} 0.35-1.00)", 0.05);
                fc->SetMinimum(0); fc->SetMaximum(1.12); fc->Draw("axis");
                const int nc = 20;
                const int gcol[6] = {kGreen + 2, kColBumpDalitz, kColEtaBkg, kColPhysBkg, kColSig, kMagenta + 2};
                for (int q = 0; q < 6; ++q) {
                    TGraph* g = new TGraph();
                    for (int k = 0; k < nc; ++k) {
                        const double c = 2.30 + 0.05 * k;
                        double v = 0;
                        if (q == 0) v = cumS(scanS[j][mEX], c);
                        else if (q == 1) v = cumS(scanS[j][mSig], c);
                        else if (q == 2) v = cumS(scanS[j][mKE], c);
                        else if (q == 3) v = cumS(scanS[j][mKM], c);
                        else if (q == 4) v = cumD(dScan[j][0], c);
                        else {
                            const double es = cumS(scanS[j][mSum], c);
                            v = es > 0 ? cumS(scanS[j][mEX], c) * scanS[j][mEX].S / (es * scanS[j][mSum].S) : 0;
                        }
                        g->SetPoint(k, c, v);
                    }
                    g->SetLineColor(gcol[q]); g->SetLineWidth(q == 0 ? 4 : 3); g->SetLineStyle(q == 5 ? 7 : 1);
                    g->SetMarkerColor(gcol[q]); g->SetMarkerStyle(q == 4 ? 20 : 1); g->SetMarkerSize(1.2);
                    g->Draw(q == 4 ? "LP same" : "L same");
                }
                VLines({gMMcut}, 0, 1.12, kRed, 2);
                TLine* l9 = new TLine(2.30, 0.9, 3.25, 0.9); l9->SetLineStyle(3); l9->Draw();
                Txt(0.15, 0.92, Form("(%s) %s: fraction kept by MM_{4} < c", j == 0 ? "c" : "f", kSel[ms[j]].tag), 0.05);
                Txt(0.17, 0.20, "dashed: sim. purity of p p #eta (X)", 0.040, kMagenta + 2);
                Txt(0.17, 0.15, "blue markers: data N - CB", 0.040, kColSig);
            }
            cv->SaveAs((fig + "/eta_enriched_mm4.png").c_str());
            cv->SaveAs((fig + "/eta_enriched_mm4.pdf").c_str());
        }
        printf("part 1 (MM4): %.1f s, cut %.2f\n", sw.RealTime(), gMMcut);
    }

    // ======================================================================
    // shared fills: M4 spectra (log 5 MeV 0.35-1.00, linear 2.5 MeV 0.40-0.70), regions
    // ======================================================================
    Binning bLog; bLog.nbin = 130; bLog.lo = 0.35; bLog.hi = 1.00;
    Binning bLin; bLin.nbin = 120; bLin.lo = 0.40; bLin.hi = 0.70;
    const Binning bReg = RegBinning();
    std::vector<std::vector<WSpec>> sLog(kNSel), sLin(kNSel), sReg(kNSel);
    std::vector<DSpec> dLog(kNSel), dLin(kNSel);
    std::vector<scan::Spectra> dReg, dPk;
    {
        TStopwatch sw;
        for (int is = 0; is < kNSel; ++is) {
            const Sel sl = SimSel(is);
            sLog[is] = FillMany(d, sl, comps, opt, bLog);
            sLin[is] = FillMany(d, sl, comps, opt, bLin);
            sReg[is] = FillMany(d, sl, comps, opt, bReg);
            dLog[is] = DataSpectrum(dd, DatSel(is), bLog.nbin, bLog.lo, bLog.hi);
            dLin[is] = DataSpectrum(dd, DatSel(is), bLin.nbin, bLin.lo, bLin.hi);
        }
        std::vector<scan::Sel> dsels;
        for (int is = 0; is < kNSel; ++is) dsels.push_back(DatSel(is));
        scan::FillOpt fo; fo.nbin = 260; fo.lo = 0.35; fo.hi = 1.00; fo.reb = 4;
        for (int k = 0; k <= kNReg; ++k) fo.windows.push_back(std::make_pair(kReg[k].lo, kReg[k].hi));
        dReg = scan::FillMany(dd, dsels, fo);          // win[2 + k] = region k
        dPk = scan::FillMany(dd, dsels, scan::FillOpt());   // standard binning for P_main
        printf("shared fills: %.1f s\n", sw.RealTime());
        // closure of the component split: sum of the six = T1 + T2
        double worst = 0;
        for (int is = 0; is < kNSel; ++is) {
            double s6 = 0; for (int c : kStack) s6 += sReg[is][c].Sall;
            worst = std::max(worst, std::fabs(s6 - sReg[is][cSum].Sall) / std::max(1e-30, sReg[is][cSum].Sall));
        }
        fprintf(fchk, "closure: max relative difference (sum of the six components - T1+T2) over the selections: %.2e\n", worst);
        // rules v2 (2026-11-13): 0.980 +- 0.041 and 5.25 +- 0.12 nb (rules v1: 0.964 +- 0.040 and 5.24 +- 0.12 nb)
        fprintf(fchk, "reference (tables/eta_engine_reference.txt, slides, peak window): peak total 0.980 +- 0.041 nb, T1+T2 5.25 +- 0.12 nb\n");
        const WSpec& r1 = sReg[sE1][cSum];
        double pt = sReg[sE1][cPC].s[3] + sReg[sE1][cPD].s[3], ptv = sReg[sE1][cPC].v[3] + sReg[sE1][cPD].v[3];
        fprintf(fchk, "this macro E1, peak window: peak total %.4f +- %.4f nb (conv. + direct, variances added), T1+T2 %.4f +- %.4f nb\n",
                pt / kLsim, std::sqrt(ptv) / kLsim, r1.s[3] / kLsim, std::sqrt(r1.v[3]) / kLsim);
    }

    // ======================================================================
    // part 2: purity, regions, M4 figures
    // ======================================================================
    std::vector<double> purAll(kNSel), purAllE(kNSel), purAllA(kNSel), purAllAE(kNSel);
    std::vector<std::vector<EffRes>> purReg(kNSel);
    // region ratios data / sim: [sel][q] default and (a)
    std::vector<std::vector<double>> rat(kNSel, std::vector<double>(kNQ)), ratE = rat, ratA = rat, ratAE = rat;
    if (parts & 2) {
        TStopwatch sw;
        // ---- purity ----
        gTab = fopen((tab + "/eta_enriched_purity.txt").c_str(), "w");
        P("# eta_enriched_purity: simulated purity of 'p p eta (X)' = (signal classes 2-5 + physical-background rows with a track from an eta)\n");
        P("#   / (T1 + T2) [default] or / (T1 + T2 + ghost-flagged true pairs) [(a)]; sum w, std + lep + eta owned rows, mode A, w (PDG)\n");
        P("# errors: event groups (pass and fail rows of one group correlated; simeng::Efficiency)\n");
        P("# M4 < 1.0: every selected row with M4 < 1.0 (the spectrum starts at 0.30)\n");
        P(Form("%-6s %-62s %16s %16s %16s %16s %16s %16s\n", "sel", "cuts (all: pre)", "M4<1.0", "M4<1.0 (a)", "below bump",
               "bump window", "peak window", "0.56-0.70"));
        for (int is = 0; is < kNSel; ++is) {
            const Sel base = SimSelM(is, -1e9, 1.0);
            const Sel pass = [](const Row& r) { return EtaX(r); };
            const EffRes e = Efficiency(d, base, pass, CompSim(), opt);
            const EffRes ea = Efficiency(d, base, pass, CompSimVar(1), opt);
            purAll[is] = e.eps; purAllE[is] = e.err; purAllA[is] = ea.eps; purAllAE[is] = ea.err;
            purReg[is] = EfficiencyVs(d, base, pass, CompSim(), opt, bReg);
            P(Form("%-6s %-62s %7.4f +- %.4f %7.4f +- %.4f %7.4f +- %.4f %7.4f +- %.4f %7.4f +- %.4f %7.4f +- %.4f\n", kSel[is].tag,
                   SelText(is).c_str(), e.eps, e.err, ea.eps, ea.err, purReg[is][0].eps, purReg[is][0].err, purReg[is][1].eps,
                   purReg[is][1].err, purReg[is][3].eps, purReg[is][3].err, purReg[is][5].eps, purReg[is][5].err));
        }
        P("\n# composition in M4 < 1.0 [nb] (default spectrum T1 + T2; ghost pairs and info classes separately)\n");
        P(Form("%-6s", "sel"));
        for (int c = 0; c < cN; ++c) if (c != cX17) P(Form(" %16s", CompShort(c)));
        P("\n");
        for (int is = 0; is < kNSel; ++is) {
            const std::vector<WSum> ws = SumMany(d, SimSelM(is, -1e9, 1.0), comps, opt);
            P(Form("%-6s", kSel[is].tag));
            for (int c = 0; c < cN; ++c) if (c != cX17) P(Form(" %7.4f+-%-7.4f", ws[c].Xsec(), ws[c].XsecErr()));
            P("\n");
        }
        fclose(gTab); gTab = nullptr;

        // ---- regions ----
        gTab = fopen((tab + "/eta_enriched_regions.txt").c_str(), "w");
        P("# eta_enriched_regions: data against simulation per M4 region (exact window sums)\n");
        P(Form("# data: rows N (event-group variance), CB = 2 sqrt(N++ N--) in 10 MeV M4 bins (ScanEngine), N - CB; nb = rows / L_data, %s\n", kLdTxt));
        P(Form("# simulation: sum w / L_sim [nb], L_sim = %.2f nb^-1, std + lep + eta owned rows, mode A, w (PDG); errors from event groups\n", kLsim));
        P("# ratio = data (N - CB) [nb] / simulation [nb]; default sim = T1 + T2; (a) = T1 + T2 + ghost-flagged true pairs\n");
        P("# component shares: fraction of T1 + T2 in the region\n");
        P("# READING (review finding 6): 0.56-0.70 and 0.70-1.00 are 80-90 %% multi-pion background (column bkgMP), which SMASH does not describe;\n");
        P("#   their ratios are listed, not used to validate the simulation; they move with the ghost treatment ((a)), the vertex slice\n");
        P("#   (eta_summary.txt section 2: central slice) and L_data (unconfirmed). The ratio near 1 below the bump is a coincidence of components\n");
        P("#   (simulated bump tail too large, eta-related background too small), not a check of the eta yield. No selection reaches eta dominance\n");
        P("#   outside the bump window (eta_enriched_purity.txt).\n");
        for (int is = 0; is < kNSel; ++is) {
            P(Form("\n## %s: %s (%s)\n", kSel[is].tag, kSel[is].name, SelText(is).c_str()));
            P(Form("%-26s %9s %8s %18s %18s %18s %16s %16s | %6s %6s %6s %6s %6s %6s\n", "region", "N", "CB", "N-CB [nb]", "sim [nb]",
                   "sim (a) [nb]", "data/sim", "data/sim (a)", "bkgMP", "bkgEta", "bumpD", "bumpC", "peakC", "peakD"));
            for (int k = 0; k <= kNReg; ++k) {
                const scan::WinSum& w = dReg[is].win[2 + k];
                const double dn = w.sig() / kLdata, de = w.sigErr() / kLdata;
                double ss, sv, sa, sav;
                if (k < kNReg) { ss = sReg[is][cSum].s[k]; sv = sReg[is][cSum].v[k]; sa = sReg[is][cSumA].s[k]; sav = sReg[is][cSumA].v[k]; }
                else { ss = sReg[is][cSum].S; sv = sReg[is][cSum].V; sa = sReg[is][cSumA].S; sav = sReg[is][cSumA].V; }
                ss /= kLsim; sa /= kLsim; sv /= kLsim * kLsim; sav /= kLsim * kLsim;
                const double r = dn / ss, re = r * std::sqrt(de * de / (dn * dn) + sv / (ss * ss));
                const double ra = dn / sa, rae = ra * std::sqrt(de * de / (dn * dn) + sav / (sa * sa));
                for (int q = 0; q < kNQ; ++q)
                    if (kQ[q] == k) { rat[is][q] = r; ratE[is][q] = re; ratA[is][q] = ra; ratAE[is][q] = rae; }
                P(Form("%-26s %9.0f %8.0f %8.4f +- %-7.4f %8.4f +- %-7.4f %8.4f +- %-7.4f %7.4f +- %-6.4f %7.4f +- %-6.4f |", kReg[k].name, w.n[0], w.cb, dn,
                       de, ss, std::sqrt(sv), sa, std::sqrt(sav), r, re, ra, rae));
                for (int c : kStack) {
                    const double v = k < kNReg ? sReg[is][c].s[k] : sReg[is][c].S;
                    P(Form(" %6.3f", v / (ss * kLsim)));
                }
                P("\n");
            }
        }
        fclose(gTab); gTab = nullptr;

        // ---- M4 figures: one per pair (without OA, with OA) ----
        struct Pair { const char* tag; int a, b; const char* title; };
        const Pair pairs[5] = {{"E12", sE2, sE1, "E2 (slides without OA) and E1 (slides)"},
                               {"E3", sE3, sE3oa, "E3 (old code) without and with OA > 2#circ"},
                               {"E4", sE4, sE4oa, "E4 (eta band) without and with OA > 2#circ"},
                               {"E12mm", sE2mm, sE1mm, "E2 and E1 with the MM_{4} cut"},
                               {"E12v", sE2v, sE1v, "E2 and E1 with both pions in the vertex fit"}};
        for (const Pair& pr : pairs) {
            TCanvas* cv = new TCanvas(Form("cM4_%s", pr.tag), "", 2700, 2300);
            cv->cd();
            MkPad(Form("pH_%s", pr.tag), 0, 0.885, 1, 1, false, 0, 0, 0, 0);
            Txt(0.01, 0.84, Form("d#sigma/dM_{#pi^{+}#pi^{-}e^{+}e^{-}}: data against the full simulation, %s", pr.title), 0.14);
            Txt(0.01, 0.66, NormLine().c_str(), 0.088);
            Txt(0.01, 0.51, Form("%s; ratio pads: data (N - CB) / simulation (full markers: T1 + T2, open: variant (a))", PreLine()), 0.088);
            SharedLegend(0.0, 0.0, 1.0, 0.45, 6, 0.085);
            const int sl[2] = {pr.a, pr.b};
            // matching y ranges of the two rows
            double yLogHi = 0, yLogLo = 1e30, yLinHi = 0;
            for (int j = 0; j < 2; ++j) {
                TH1D* hN = dLog[sl[j]].Hist("us", Form("tmpN%d", j), 2);
                yLogHi = std::max(yLogHi, HMax(hN));
                yLogLo = std::min(yLogLo, std::max(HMinPos(dLog[sl[j]].Hist("cb", Form("tmpC%d", j), 2)) * 0.3, yLogHi * 1e-4));
                TH1D* hS = dLin[sl[j]].Hist("sig", Form("tmpS%d", j), 2);
                TH1D* hM = sLin[sl[j]][cSumA].Hist(Form("tmpM%d", j), 2);
                yLinHi = std::max(yLinHi, std::max(HMax(hS), HMax(hM)));
            }
            for (int j = 0; j < 2; ++j) {
                const int is = sl[j];
                const double y1 = j == 0 ? 0.4425 : 0.0, y2 = j == 0 ? 0.885 : 0.4425;
                const double ysplit = y1 + 0.30 * (y2 - y1);
                for (int v = 0; v < 2; ++v) {
                    const bool lin = v == 1;
                    const double x1 = lin ? 0.5 : 0.0, x2 = lin ? 1.0 : 0.5;
                    const std::vector<WSpec>& S = lin ? sLin[is] : sLog[is];
                    const DSpec& D = lin ? dLin[is] : dLog[is];
                    const std::string ids = Form("%s_%d_%d", pr.tag, j, v);
                    const char* id = ids.c_str();
                    cv->cd();
                    MkPad(Form("pm_%s", id), x1, ysplit, x2, y2, !lin, 0.12, 0.02, 0.11, 0.0);
                    TH1D* hN = D.Hist("us", Form("hN_%s", id), 2);
                    TH1D* hC = D.Hist("cb", Form("hC_%s", id), 2);
                    TH1D* hS = D.Hist("sig", Form("hS_%s", id), 2);
                    TH1D* hSum = S[cSum].Hist(Form("hSum_%s", id), 2);
                    TH1D* hSumA = S[cSumA].Hist(Form("hSumA_%s", id), 2);
                    TH1D* fr = (TH1D*)hSum->Clone(Form("fr_%s", id)); fr->Reset(); Axes(fr, kXM4, kYM4, 0.055, 1.05);
                    if (lin) { fr->SetMinimum(0); fr->SetMaximum(yLinHi * 1.15); }
                    else { fr->SetMinimum(yLogLo); fr->SetMaximum(yLogHi * 60); }
                    fr->Draw("axis");
                    if (lin) {
                        TH1D* hs[6];
                        for (int k = 0; k < 6; ++k) {
                            hs[k] = S[kStack[k]].Hist(Form("hst_%s_%d", id, k), 2);
                            if (k > 0) hs[k]->Add(hs[k - 1]);
                            StyleH(hs[k], kStackLine[k], 1, 1, kStackFill[k]);
                        }
                        for (int k = 5; k >= 0; --k) hs[k]->Draw("hist same");
                        TH1D* hx = S[cX17].Hist(Form("hx_%s", id), 2); StyleH(hx, kColX17, 3, 2); hx->Draw("hist same");
                    } else {
                        TH1D* hk = S[cKM].Hist(Form("hkm_%s", id), 2); StyleH(hk, kColPhysBkg, 1, 1, kColPhysFill); hk->Draw("hist same");
                        for (int k = 1; k < 6; ++k) {
                            TH1D* h = S[kStack[k]].Hist(Form("hl_%s_%d", id, k), 2); StyleH(h, kStackLine[k], 3, 1); h->Draw("hist same");
                        }
                        StyleM(hN, kColData, 20, 0.9); hN->Draw("E1 same");
                        StyleM(hC, kColCB, 21, 0.9); hC->Draw("E1 same");
                    }
                    StyleH(hSum, kColSim, 4, 1); hSum->Draw("hist same");
                    StyleH(hSumA, kColSim, 3, 2); hSumA->Draw("hist same");
                    StyleM(hS, kColSig, 20, lin ? 1.1 : 0.9); hS->Draw("E1 same");
                    VLines({scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}, fr->GetMinimum(), fr->GetMaximum());
                    if (!lin) VLines({0.56, 0.70}, fr->GetMinimum(), fr->GetMaximum(), kGray + 2, 2);
                    fr->Draw("axis same");
                    Txt(0.13, 0.91, Form("(%s) %s; %s", j == 0 ? (lin ? "b" : "a") : (lin ? "d" : "c"), SelLatex(is).c_str(),
                                         lin ? "linear, 2.5 MeV, stacked" : "log, 5 MeV"), 0.040);
                    if (!lin) {
                        Txt(0.15, 0.83, "data/sim:", 0.040);
                        const char* qn[kNQ] = {"0.35-0.455", "bump", "peak", "0.56-0.70", "0.70-1.00", "0.35-1.00"};
                        for (int q = 0; q < kNQ; ++q)
                            Txt(0.28 + 0.24 * (q % 3), 0.83 - 0.055 * (q / 3), Form("%s: %.2f#pm%.2f, (a) %.2f", qn[q], rat[is][q], ratE[is][q], ratA[is][q]), 0.030);
                        Txt(0.15, 0.725, Form("purity p p #eta (X), M_{4} < 1: %.3f #pm %.3f", purAll[is], purAllE[is]), 0.036, kColPeakTotal);
                    } else {
                        Txt(0.15, 0.83, Form("sim. peak total (window): %.3f nb; X17 (BR 10^{-5}) dashed red, not stacked",
                                             (S[cPC].Window(scan::kPeakLo, scan::kPeakHi).first + S[cPD].Window(scan::kPeakLo, scan::kPeakHi).first) / kLsim), 0.034);
                    }
                    // ratio pad
                    cv->cd();
                    MkPad(Form("pr_%s", id), x1, y1, x2, ysplit, false, 0.12, 0.02, 0.0, 0.30);
                    TH1D* r1 = RatioH(hS, hSum, Form("r1_%s", id));
                    TH1D* r2 = RatioH(hS, hSumA, Form("r2_%s", id));
                    TH1D* frr = (TH1D*)r1->Clone(Form("frr_%s", id)); frr->Reset();
                    Axes(frr, kXM4, "data/sim", 0.12, 0.45);
                    frr->GetYaxis()->SetNdivisions(505);
                    frr->SetMinimum(0.0); frr->SetMaximum(2.6); frr->Draw("axis");
                    TLine* l1 = new TLine(frr->GetXaxis()->GetXmin(), 1, frr->GetXaxis()->GetXmax(), 1); l1->SetLineStyle(2); l1->Draw();
                    StyleM(r1, kColSig, 20, 0.8); r1->Draw("E1 same");
                    StyleM(r2, kColSim, 24, 0.8); r2->Draw("E1 same");
                    VLines({scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}, 0, 2.6);
                }
            }
            cv->SaveAs((fig + "/eta_enriched_m4_" + pr.tag + ".png").c_str());
            cv->SaveAs((fig + "/eta_enriched_m4_" + pr.tag + ".pdf").c_str());
        }

        // ---- summary figure: purity and region ratios per selection ----
        {
            TCanvas* cv = new TCanvas("cSumm", "", 2700, 1500);
            cv->cd();
            MkPad("pHs", 0, 0.88, 1, 1, false, 0, 0, 0, 0);
            Txt(0.01, 0.75, "Eta-enriched selections: simulated purity of p p #eta (X) and data/simulation per M_{4} region", 0.20);
            Txt(0.01, 0.45, NormLine().c_str(), 0.12);
            Txt(0.01, 0.17, Form("%s; selections defined in tables/eta_enriched_regions.txt (MM_{4} cut %.2f GeV/c^{2})", PreLine(), gMMcut), 0.12);
            for (int p = 0; p < 3; ++p) {
                cv->cd();
                MkPad(Form("pS%d", p), p / 3., 0.0, (p + 1) / 3., 0.88, false, 0.13, 0.03, 0.09, 0.14);
                TH1D* fr = new TH1D(Form("sfr%d", p), "", kNSel, 0, kNSel); fr->SetDirectory(nullptr);
                for (int is = 0; is < kNSel; ++is) fr->GetXaxis()->SetBinLabel(is + 1, kSel[is].tag);
                Axes(fr, "", p == 0 ? "simulated purity of p p #eta (X)" : "data (N - CB) / simulation", 0.05, 1.2);
                fr->GetXaxis()->SetLabelSize(0.050); fr->GetXaxis()->LabelsOption("h");
                fr->SetMinimum(p == 0 ? 0 : 0.5); fr->SetMaximum(p == 0 ? 1.05 : 2.0); fr->Draw("axis");
                if (p > 0) { TLine* l = new TLine(0, 1, kNSel, 1); l->SetLineStyle(2); l->Draw(); }
                TLegend* lg = new TLegend(0.16, p == 0 ? 0.70 : 0.66, 0.97, 0.90);
                lg->SetNColumns(2); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.035);
                if (p == 0) {
                    const int col[5] = {kBlack, kColSim, kColBumpDalitz, kColPeakDirect, kColPhysBkg};
                    const char* nm[5] = {"M_{4} < 1.0", "M_{4} < 1.0, variant (a)", "bump window", "peak window", "0.56-0.70"};
                    for (int q = 0; q < 5; ++q) {
                        TGraphErrors* g = new TGraphErrors();
                        for (int is = 0; is < kNSel; ++is) {
                            double v, e;
                            if (q == 0) { v = purAll[is]; e = purAllE[is]; }
                            else if (q == 1) { v = purAllA[is]; e = purAllAE[is]; }
                            else { const int k = q == 2 ? 1 : q == 3 ? 3 : 5; v = purReg[is][k].eps; e = purReg[is][k].err; }
                            g->SetPoint(is, is + 0.2 + 0.15 * q, v); g->SetPointError(is, 0, e);
                        }
                        g->SetMarkerStyle(q == 1 ? 24 : 20 + (q > 1 ? q : 0)); g->SetMarkerColor(col[q]); g->SetLineColor(col[q]);
                        g->SetMarkerSize(1.5); g->SetLineWidth(2); g->Draw("P same");
                        lg->AddEntry(g, nm[q], "p");
                    }
                    Txt(0.15, 0.93, "(a) simulated purity of p p #eta (X)", 0.042);
                    Txt(0.16, 0.13, "purity = (signal classes 2-5 + eta-related bkg) / (T1 + T2)", 0.032);
                } else {
                    const int col[kNQ] = {kOrange + 7, kColBumpDalitz, kColPeakDirect, kColPhysBkg, kGray + 3, kBlack};
                    const char* nm[kNQ] = {"0.35-0.455", "bump window", "peak window", "0.56-0.70", "0.70-1.00", "0.35-1.00"};
                    for (int q = 0; q < kNQ; ++q) {
                        TGraphErrors* g = new TGraphErrors();
                        for (int is = 0; is < kNSel; ++is) {
                            g->SetPoint(is, is + 0.15 + 0.13 * q, p == 1 ? rat[is][q] : ratA[is][q]);
                            g->SetPointError(is, 0, p == 1 ? ratE[is][q] : ratAE[is][q]);
                        }
                        g->SetMarkerStyle(20 + q); g->SetMarkerColor(col[q]); g->SetLineColor(col[q]); g->SetMarkerSize(1.5);
                        g->SetLineWidth(2); g->Draw("P same");
                        lg->AddEntry(g, nm[q], "p");
                    }
                    Txt(0.15, 0.93, p == 1 ? "(b) data/sim, simulation T1 + T2" : "(c) data/sim, variant (a): T1 + T2 + ghost pairs", 0.042);
                }
                lg->Draw();
            }
            cv->SaveAs((fig + "/eta_enriched_summary.png").c_str());
            cv->SaveAs((fig + "/eta_enriched_summary.pdf").c_str());
        }
        printf("part 2: %.1f s\n", sw.RealTime());
    }

    // ======================================================================
    // part 4: template fits
    // ======================================================================
    if (parts & 4) {
        TStopwatch sw;
        const std::vector<Comp> fc = FitComps();
        struct Rng { double lo, hi; int nb; };
        const Rng rng[2] = {{0.40, 0.75, 70}, {0.43, 0.70, 54}};
        // models: groups of fit components, per ghost variant (0 default, 1 (a))
        struct Model { const char* tag; const char* desc; std::vector<int> g[2]; std::vector<std::string> pn; };
        std::vector<Model> models = {
            {"M1", "bump, peak (direct+conv fixed), bkg", {{fBump, fPeak, fK}, {fBump, fPeak, fKA}}, {"a_bump", "a_peak", "a_bkg"}},
            {"M2", "bump, peak direct, peak conv, bkg", {{fBump, fPD, fPC, fK}, {fBump, fPD, fPC, fKA}}, {"a_bump", "a_pdir", "a_pconv", "a_bkg"}},
            {"M3", "bump, peak, bkg eta-rel, bkg multi-pion", {{fBump, fPeak, fKE, fKM}, {fBump, fPeak, fKEA, fKMA}}, {"a_bump", "a_peak", "a_bkgEta", "a_bkgMP"}},
            {"M4", "bump, peak direct, peak conv, bkg eta-rel, bkg multi-pion", {{fBump, fPD, fPC, fKE, fKM}, {fBump, fPD, fPC, fKEA, fKMA}},
             {"a_bump", "a_pdir", "a_pconv", "a_bkgEta", "a_bkgMP"}}};
        // simulated contents (all M4 < 1.1, and windows) in expected data rows
        const double toRows = kLdata / kLsim;
        // P_main on the data, standard ScanEngine binning, vz > -160 (this selection)
        std::vector<scan::PeakFit> pm(kNSel);
        for (int is = 0; is < kNSel; ++is) pm[is] = scan::FitPeakMain(dPk[is]);
        // P_main applied to the simulation (pseudo-data = simulated spectrum in expected data rows, template variances)
        //   0: T1+T2, sim M4 x 0.9952, sigma 5.4 MeV (data shape)   1: the same with sigma 6.2 MeV (simulated width)
        //   2: variant (a), x 0.9952, 5.4 MeV                       3: peak template alone (direct + conv.), x 0.9952, 5.4 MeV
        //   4: T1+T2 unscaled, Gaussian at the simulated position 0.5460, 6.2 MeV
        // pmD: the same hand-made fit on the data (diagonal N + CB variance), cross-check of scan::FitPeakMain
        const int kNPM = 5;
        std::vector<std::vector<PM>> pmS(kNSel, std::vector<PM>(kNPM));
        std::vector<PM> pmD(kNSel);
        for (int is = 0; is < kNSel; ++is) {
            const scan::Spectra& sp = dPk[is];
            std::vector<double> vd(sp.nbin);
            for (int i = 0; i < sp.nbin; ++i) vd[i] = sp.n[0][i] + sp.cbVar[i];
            pmD[is] = PmainLin(sp.sig, vd, sp.lo, sp.bw, scan::kPeakMean, scan::kPeakSigma);
            for (int sc = 0; sc < 2; ++sc) {
                Binning bP; bP.nbin = 240; bP.lo = 0.35; bP.hi = 0.95;
                const double scl = sc ? kScaleK0S : 1.;
                bP.x = [scl](const Row& r) { return r.m4 * scl; };
                const std::vector<WSpec> T = FillMany(d, SimSel(is), {CompSim(), CompSimVar(1), CompPeakTotal()}, opt, bP);
                std::vector<std::vector<double>> y(3, std::vector<double>(bP.nbin)), v = y;
                for (int c = 0; c < 3; ++c)
                    for (int i = 0; i < bP.nbin; ++i) { y[c][i] = T[c].s[i] * toRows; v[c][i] = T[c].v[i] * toRows * toRows; }
                if (sc) {
                    pmS[is][0] = PmainLin(y[0], v[0], bP.lo, bP.bw(), scan::kPeakMean, scan::kPeakSigma);
                    pmS[is][1] = PmainLin(y[0], v[0], bP.lo, bP.bw(), scan::kPeakMean, 0.0062);
                    pmS[is][2] = PmainLin(y[1], v[1], bP.lo, bP.bw(), scan::kPeakMean, scan::kPeakSigma);
                    pmS[is][3] = PmainLin(y[2], v[2], bP.lo, bP.bw(), scan::kPeakMean, scan::kPeakSigma);
                } else
                    pmS[is][4] = PmainLin(y[0], v[0], bP.lo, bP.bw(), 0.5460, 0.0062);
            }
        }
        auto pmRatio = [&](int is, int k, double& r, double& e) {
            r = pm[is].S / pmS[is][k].S;
            e = std::fabs(r) * std::sqrt(std::pow(pm[is].eS / pm[is].S, 2) + std::pow(pmS[is][k].eS / pmS[is][k].S, 2));
        };
        // data-study numbers (reports2/reconcile.md; P_main, NO vz cut): slides 3168 +- 304 +- 717, old code 3650 +- 314 +- 700
        gTab = fopen((tab + "/eta_enriched_fit.txt").c_str(), "w");
        P("# eta_enriched_fit: template fits of data N - CB with free factors on the SMASH templates (a comparison for discussion, NOT a measurement)\n");
        P(Form("# data: N - CB in 5 MeV bins, covariance = Poisson N + interpolated CB (10 MeV); templates in expected data rows (x %s / L_sim)\n", kLdTxt));
        P("# chi2 = (y - sum a_k T_k)^T V^-1 (...), V = data covariance + diag(sum a_k^2 Var T_k), iterated; template bin-to-bin group correlations neglected\n");
        P("# a = 1 means the SMASH prediction with L_data = 6080.5 nb^-1 (unconfirmed): every a scales with 6080.5 / L_data\n");
        P("# models: M1 bump (Dalitz + conv., SMASH ratio), peak (direct + conv., SMASH ratio), physical bkg; M2 peak direct and conv. free;\n");
        P("#   M3 bkg eta-related and multi-pion free; M4 both splits. variant (a): bkg + ghost-flagged true pairs. scale: sim M4 x 0.9952 (K0S)\n");
        P("# implied data contents: peak = a x (simulated peak, all M4) in rows; bump = a x (simulated bump, all M4) in rows\n");
        P("# P_main (this macro): ScanEngine FitPeakMain on the data of the same selection (vz > -160): Gaussian 0.5439 / 5.4 MeV frozen, pol2, 0.5225-0.61\n");
        P("# data study (reports2/reconcile.md, P_main, NO vz cut): slides 3168 +- 304 (stat) +- 717 (model); old code 3650 +- 314 +- 700; no cut 4657 +- 390 +- 1233\n");
        P("# REVIEW (finding 5): the fits are poor (chi2/ndf up to 3.3, p down to 1e-14); the errors below are NOT scaled (multiply by sqrt(chi2/ndf));\n");
        P("#   the model spread (M1/M3, range, scale, ghost variant) is far larger than the errors. The factors compare shapes; they are NOT yields:\n");
        P("#   alpha_bump lies between about 0.12 and 0.50 over all M1/M3 variants (not 'stable'); a_bkgEta > 1 together with alpha_bump < 1 is\n");
        P("#   contradictory. Never combine the bump and peak factors into an eta production factor. 'implied' contents are not measurements.\n");
        P("# REVIEW (finding 1): the peak data/sim depends on the method (P_main on both sides, template fit, side bands, all-M4 content):\n");
        P("#   see tables/eta_summary_peak_methods.txt\n");
        P("\n## peak by the data method: P_main on the data and on the simulation (expected data rows; sim errors from the template variances)\n");
        P("# sim M4 x 0.9952 (K0S scale) so that the simulated peak sits at the data position; Gaussian 0.5439, sigma 5.4 MeV (data) or 6.2 MeV (sim)\n");
        P("# 'unscaled': sim M4 as produced, Gaussian at the simulated position 0.5460, 6.2 MeV; 'peak only': the fit on the peak template alone\n");
        P(Form("%-6s %16s %16s | %16s %16s %16s %16s %16s | %14s %14s %14s %14s\n", "sel", "data P_main", "data (hand)", "sim 5.4", "sim 6.2",
               "sim (a) 5.4", "peak only 5.4", "sim unscaled", "data/sim 5.4", "data/sim 6.2", "d/s (a) 5.4", "d/s unscaled"));
        for (int is = 0; is < kNSel; ++is) {
            double r[4], e[4];
            pmRatio(is, 0, r[0], e[0]); pmRatio(is, 1, r[1], e[1]); pmRatio(is, 2, r[2], e[2]); pmRatio(is, 4, r[3], e[3]);
            P(Form("%-6s %7.0f +- %-6.0f %7.0f +- %-6.0f | %7.0f +- %-6.0f %7.0f +- %-6.0f %7.0f +- %-6.0f %7.0f +- %-6.0f %7.0f +- %-6.0f | %5.3f +- %-5.3f %5.3f +- %-5.3f %5.3f +- %-5.3f %5.3f +- %-5.3f\n",
                   kSel[is].tag, pm[is].S, pm[is].eS, pmD[is].S, pmD[is].eS, pmS[is][0].S, pmS[is][0].eS, pmS[is][1].S, pmS[is][1].eS,
                   pmS[is][2].S, pmS[is][2].eS, pmS[is][3].S, pmS[is][3].eS, pmS[is][4].S, pmS[is][4].eS, r[0], e[0], r[1], e[1], r[2],
                   e[2], r[3], e[3]));
        }
        P(Form("chi2/ndf of the simulated P_main fits (5.4 MeV, scaled): "));
        for (int is = 0; is < kNSel; ++is) P(Form("%s %.0f/%d  ", kSel[is].tag, pmS[is][0].chi2, pmS[is][0].ndf));
        P("\n");
        struct Store { FitRes r; double pkRows, pkRowsE, bmRows, bmRowsE; };
        // store[is][range][scale][ghost][model]
        std::vector<std::vector<std::vector<std::vector<std::vector<Store>>>>> st(
            kNSel, std::vector<std::vector<std::vector<std::vector<Store>>>>(
                       2, std::vector<std::vector<std::vector<Store>>>(2, std::vector<std::vector<Store>>(2, std::vector<Store>(4)))));
        // keep the spectra of E2 and E1 (range 0, no scale) for the figure
        std::vector<WSpec> keepT[2]; DSpec keepD[2];
        for (int is = 0; is < kNSel; ++is) {
            const std::vector<WSum> tot = SumMany(d, SimSel(is), fc, opt);
            P(Form("\n## %s: %s (%s)\n", kSel[is].tag, kSel[is].name, SelText(is).c_str()));
            P(Form("simulated contents (all M4, expected data rows): bump %.0f +- %.0f, peak total %.0f +- %.0f (direct %.0f, conv. %.0f)\n",
                   tot[fBump].S * toRows, std::sqrt(tot[fBump].V) * toRows, tot[fPeak].S * toRows, std::sqrt(tot[fPeak].V) * toRows,
                   tot[fPD].S * toRows, tot[fPC].S * toRows));
            P(Form("data P_main (vz > -160, this selection): S = %.0f +- %.0f (stat), chi2/ndf %.1f/%d\n", pm[is].S, pm[is].eS, pm[is].chi2,
                   pm[is].ndf));
            P(Form("%-3s %-10s %-5s %-3s %-58s %8s %6s %22s %22s\n", "rng", "range", "scale", "var", "factors a (+- error)", "chi2/ndf", "p",
                   "implied peak [rows]", "implied bump [rows]"));
            for (int ir = 0; ir < 2; ++ir) {
                const DSpec D = DataSpectrum(dd, DatSel(is), rng[ir].nb, rng[ir].lo, rng[ir].hi);
                std::vector<int> bins;
                for (int i = 0; i < D.nbin; ++i) bins.push_back(i);
                for (int sc = 0; sc < 2; ++sc) {
                    Binning b; b.nbin = rng[ir].nb; b.lo = rng[ir].lo; b.hi = rng[ir].hi;
                    const double scl = sc ? kScaleK0S : 1.;
                    b.x = [scl](const Row& r) { return r.m4 * scl; };
                    const std::vector<WSpec> T = FillMany(d, SimSel(is), fc, opt, b);
                    if (ir == 0 && sc == 0 && (is == sE1 || is == sE2)) { keepT[is == sE1 ? 1 : 0] = T; keepD[is == sE1 ? 1 : 0] = D; }
                    for (int gv = 0; gv < 2; ++gv)
                        for (int im = 0; im < 4; ++im) {
                            const std::vector<int>& grp = models[im].g[gv];
                            std::vector<std::vector<double>> TT, TV;
                            for (int k : grp) {
                                std::vector<double> t(b.nbin), tv(b.nbin);
                                for (int i = 0; i < b.nbin; ++i) { t[i] = T[k].s[i] * toRows; tv[i] = T[k].v[i] * toRows * toRows; }
                                TT.push_back(t); TV.push_back(tv);
                            }
                            Store so;
                            so.r = TemplateFit(D.sp, bins, TT, TV);
                            // implied contents
                            std::vector<double> cp(grp.size(), 0.), cbm(grp.size(), 0.);
                            for (size_t k = 0; k < grp.size(); ++k) {
                                if (grp[k] == fPeak || grp[k] == fPD || grp[k] == fPC) cp[k] = tot[grp[k]].S * toRows;
                                if (grp[k] == fBump) cbm[k] = tot[grp[k]].S * toRows;
                            }
                            double pk = 0, pkv = 0, bm = 0, bmv = 0;
                            for (size_t k = 0; k < grp.size(); ++k) {
                                pk += so.r.a[k] * cp[k]; bm += so.r.a[k] * cbm[k];
                                for (size_t l = 0; l < grp.size(); ++l) { pkv += cp[k] * cp[l] * so.r.cov[k][l]; bmv += cbm[k] * cbm[l] * so.r.cov[k][l]; }
                            }
                            so.pkRows = pk; so.pkRowsE = std::sqrt(std::max(0., pkv)); so.bmRows = bm; so.bmRowsE = std::sqrt(std::max(0., bmv));
                            st[is][ir][sc][gv][im] = so;
                            std::string pars;
                            for (size_t k = 0; k < grp.size(); ++k)
                                pars += Form("%s %.3f+-%.3f ", models[im].pn[k].c_str(), so.r.a[k], so.r.ea[k]);
                            P(Form("%-3d %.2f-%.2f %-5s %-3s %s %-58s %5.1f/%-3d %6.3g %10.0f +- %-9.0f %10.0f +- %-9.0f\n", ir, rng[ir].lo, rng[ir].hi,
                                   sc ? "K0S" : "none", gv ? "(a)" : "T12", models[im].tag, pars.c_str(), so.r.chi2, so.r.ndf,
                                   TMath::Prob(so.r.chi2, so.r.ndf), so.pkRows, so.pkRowsE, so.bmRows, so.bmRowsE));
                        }
                }
            }
        }
        fclose(gTab); gTab = nullptr;

        // ---- figure: fits of E2 and E1 (range 0.40-0.75, no scale, default variant), M1 and M4 ----
        {
            TCanvas* cv = new TCanvas("cFit", "", 2700, 2200);
            cv->cd();
            MkPad("pHf", 0, 0.89, 1, 1, false, 0, 0, 0, 0);
            Txt(0.01, 0.82, "Template fit of data N - CB with free factors on the SMASH components, M_{4} 0.40-0.75 GeV/c^{2} (comparison for discussion, not a measurement)", 0.14);
            Txt(0.01, 0.64, NormLine().c_str(), 0.09);
            Txt(0.01, 0.48, "#chi^{2}: data covariance (N + interpolated CB) + #alpha^{2} #times template variance; components drawn #times their fitted #alpha; pulls = (data - fit) / #sigma_{diag}", 0.09);
            SharedLegend(0.0, 0.0, 1.0, 0.42, 5, 0.085, true, false, false);
            const int im2[2] = {0, 2};
            const int sl[2] = {sE2, sE1};
            double ymax = 0;
            for (int j = 0; j < 2; ++j) ymax = std::max(ymax, HMax(keepD[j].Hist("sig", Form("tmpF%d", j), 2)));
            for (int j = 0; j < 2; ++j) {
                const int is = sl[j];
                // components of the figure at the binning of the fit (5 MeV, 0.40-0.75)
                Binning b; b.nbin = rng[0].nb; b.lo = rng[0].lo; b.hi = rng[0].hi;
                const std::vector<WSpec> C = FillMany(d, SimSel(is), comps, opt, b);
                for (int m = 0; m < 2; ++m) {
                    const int im = im2[m];
                    const Store& so = st[is][0][0][0][im];
                    const double x1 = m * 0.5, x2 = x1 + 0.5, y1 = j == 0 ? 0.445 : 0.0, y2 = j == 0 ? 0.89 : 0.445;
                    const double ys = y1 + 0.27 * (y2 - y1);
                    const std::string ids = Form("f%d%d", j, m);
                    const char* id = ids.c_str();
                    cv->cd();
                    MkPad(Form("pf_%s", id), x1, ys, x2, y2, false, 0.12, 0.02, 0.11, 0.0);
                    // alpha per drawn component
                    double al[cN];
                    for (int c = 0; c < cN; ++c) al[c] = 0;
                    if (im == 0) { al[cBD] = al[cBC] = so.r.a[0]; al[cPC] = al[cPD] = so.r.a[1]; al[cKE] = al[cKM] = so.r.a[2]; }
                    else if (im == 2) { al[cBD] = al[cBC] = so.r.a[0]; al[cPC] = al[cPD] = so.r.a[1]; al[cKE] = so.r.a[2]; al[cKM] = so.r.a[3]; }
                    else { al[cBD] = al[cBC] = so.r.a[0]; al[cPD] = so.r.a[1]; al[cPC] = so.r.a[2]; al[cKE] = so.r.a[3]; al[cKM] = so.r.a[4]; }
                    TH1D* hs[6];
                    for (int k = 0; k < 6; ++k) {
                        hs[k] = C[kStack[k]].Hist(Form("fs_%s_%d", id, k), 2);
                        hs[k]->Scale(al[kStack[k]]);
                        if (k > 0) hs[k]->Add(hs[k - 1]);
                        StyleH(hs[k], kStackLine[k], 1, 1, kStackFill[k]);
                    }
                    TH1D* hD = keepD[j].Hist("sig", Form("fD_%s", id), 2);
                    TH1D* fr = (TH1D*)hD->Clone(Form("ffr_%s", id)); fr->Reset(); Axes(fr, kXM4, kYM4, 0.055, 1.05);
                    fr->SetMinimum(0); fr->SetMaximum(ymax * 1.18); fr->Draw("axis");
                    for (int k = 5; k >= 0; --k) hs[k]->Draw("hist same");
                    TH1D* hSum = (TH1D*)hs[5]->Clone(Form("fsum_%s", id)); StyleH(hSum, kColSim, 4, 1); hSum->SetFillStyle(0); hSum->Draw("hist same");
                    StyleM(hD, kColSig, 20, 1.1); hD->Draw("E1 same");
                    VLines({scan::kBumpLo, scan::kBumpHi, scan::kPeakLo, scan::kPeakHi}, 0, ymax * 1.18);
                    fr->Draw("axis same");
                    Txt(0.13, 0.91, Form("(%s) %s; model %s", j == 0 ? (m ? "b" : "a") : (m ? "d" : "c"), SelLatex(is).c_str(), models[im].tag), 0.036);
                    std::string pars;
                    for (int k = 0; k < so.r.np; ++k) pars += Form("%s = %.2f #pm %.2f   ", models[im].pn[k].c_str(), so.r.a[k], so.r.ea[k]);
                    Txt(0.15, 0.84, pars.c_str(), 0.034);
                    Txt(0.15, 0.785, Form("#chi^{2}/ndf = %.1f/%d; implied peak %.0f #pm %.0f rows (data P_{main}: %.0f #pm %.0f)", so.r.chi2, so.r.ndf, so.pkRows,
                                          so.pkRowsE, pm[is].S, pm[is].eS), 0.034);
                    // pulls
                    cv->cd();
                    MkPad(Form("pp_%s", id), x1, y1, x2, ys, false, 0.12, 0.02, 0.0, 0.32);
                    TH1D* pu = (TH1D*)hD->Clone(Form("pu_%s", id)); pu->Reset();
                    for (int i = 1; i <= pu->GetNbinsX(); ++i) {
                        double tv = 0;
                        for (int k = 0; k < 6; ++k) {
                            const double e = C[kStack[k]].v[i - 1] > 0 ? std::sqrt(C[kStack[k]].v[i - 1]) / (kLsim * b.bw()) * al[kStack[k]] : 0;
                            tv += e * e;
                        }
                        const double den = std::sqrt(hD->GetBinError(i) * hD->GetBinError(i) + tv);
                        if (den > 0) pu->SetBinContent(i, (hD->GetBinContent(i) - hs[5]->GetBinContent(i)) / den);
                    }
                    Axes(pu, kXM4, "pull", 0.12, 0.45);
                    pu->GetYaxis()->SetNdivisions(505);
                    pu->SetMinimum(-4.5); pu->SetMaximum(4.5);
                    pu->SetFillColor(kColSig); pu->SetLineColor(kColSig); pu->SetFillStyle(1001);
                    pu->Draw("hist");
                    TLine* l0 = new TLine(b.lo, 0, b.hi, 0); l0->Draw();
                }
            }
            cv->SaveAs((fig + "/eta_enriched_fit_E12.png").c_str());
            cv->SaveAs((fig + "/eta_enriched_fit_E12.pdf").c_str());
        }
        // ---- figure: alpha values across selections and variants ----
        {
            TCanvas* cv = new TCanvas("cAl", "", 2700, 1900);
            cv->cd();
            MkPad("pHa", 0, 0.90, 1, 1, false, 0, 0, 0, 0);
            Txt(0.01, 0.72, "Fitted factors #alpha (data / SMASH) per selection: stability across fit ranges, ghost variant and mass scale", 0.20);
            Txt(0.01, 0.40, Form("#alpha = 1: SMASH with %s; every #alpha scales with 6080.5 / L_{data}. Fit range 0.40-0.75 unless stated; errors: data + template variance", kLdTex), 0.13);
            Txt(0.01, 0.10, "markers: black full = M1, 0.40-0.75; black open = 0.43-0.70; teal square = variant (a); orange triangle = sim M_{4} #times 0.9952 (K0S scale)", 0.13);
            struct Pan { const char* t; int im; int k; int k2; };
            const Pan pans[6] = {{"#alpha_{bump} (M1)", 0, 0, -1}, {"#alpha_{peak} (M1)", 0, 1, -1}, {"#alpha_{bkg} (M1)", 0, 2, -1},
                                 {"#chi^{2}/ndf (M1)", 0, -1, -1}, {"peak factor from P_{main}: data / simulation", 0, 1, -1},
                                 {"M3: #alpha_{bkg eta-rel.} (full), #alpha_{bkg multi-pion} (open)", 2, 2, 3}};
            for (int p = 0; p < 6; ++p) {
                cv->cd();
                const double x1 = (p % 3) / 3., x2 = x1 + 1. / 3., y1 = p < 3 ? 0.45 : 0.0, y2 = p < 3 ? 0.90 : 0.45;
                MkPad(Form("pa%d", p), x1, y1, x2, y2, false, 0.13, 0.03, 0.10, 0.11);
                double lo = 1e30, hi = -1e30;
                auto upd = [&](double v, double e) { lo = std::min(lo, v - e); hi = std::max(hi, v + e); };
                for (int is = 0; is < kNSel; ++is)
                    for (int ir = 0; ir < 2; ++ir)
                        for (int sc = 0; sc < 2; ++sc)
                            for (int gv = 0; gv < 2; ++gv) {
                                const Store& so = st[is][ir][sc][gv][pans[p].im];
                                if (!so.r.ok) continue;
                                if (pans[p].k < 0) upd(so.r.chi2 / so.r.ndf, 0);
                                else {
                                    upd(so.r.a[pans[p].k], so.r.ea[pans[p].k]);
                                    if (pans[p].k2 >= 0) upd(so.r.a[pans[p].k2], so.r.ea[pans[p].k2]);
                                }
                            }
                lo = std::max(std::min(lo, 0.), -1.0); hi = std::min(std::max(hi, 1.2), 4.0);
                if (p == 4) { lo = 0; hi = 1.6; }
                TH1D* fr = new TH1D(Form("afr%d", p), "", kNSel, 0, kNSel); fr->SetDirectory(nullptr);
                for (int is = 0; is < kNSel; ++is) fr->GetXaxis()->SetBinLabel(is + 1, kSel[is].tag);
                Axes(fr, "", pans[p].k < 0 ? "#chi^{2}/ndf" : "#alpha", 0.05, 1.2);
                fr->GetXaxis()->SetLabelSize(0.050);
                fr->GetXaxis()->LabelsOption("h");
                fr->SetMinimum(std::min(0., lo - 0.1 * (hi - lo))); fr->SetMaximum(hi + 0.25 * (hi - lo)); fr->Draw("axis");
                TLine* l = new TLine(0, 1, kNSel, 1); l->SetLineStyle(2); l->Draw();
                if (p == 4) {
                    // P_main ratios: full 5.4 MeV, open 6.2 MeV, teal (a), grey unscaled; orange: alpha_peak of M1 with the K0S scale
                    const int ks[4] = {0, 1, 2, 4};
                    const int mk[4] = {20, 24, 21, 25}, cl[4] = {kBlack, kBlack, kColSim, kGray + 2};
                    for (int q = 0; q < 5; ++q) {
                        TGraphErrors* g = new TGraphErrors();
                        for (int is = 0; is < kNSel; ++is) {
                            double r, e;
                            if (q < 4) pmRatio(is, ks[q], r, e);
                            else { r = st[is][0][1][0][0].r.a[1]; e = st[is][0][1][0][0].r.ea[1]; }
                            g->SetPoint(is, is + 0.15 + 0.17 * q, r); g->SetPointError(is, 0, e);
                        }
                        g->SetMarkerStyle(q < 4 ? mk[q] : 22); g->SetMarkerColor(q < 4 ? cl[q] : kOrange + 7);
                        g->SetLineColor(q < 4 ? cl[q] : kOrange + 7); g->SetMarkerSize(1.4); g->SetLineWidth(2); g->Draw("P same");
                    }
                    Txt(0.15, 0.92, "(e) peak: P_{main}(data) / P_{main}(sim, M_{4}#times0.9952)", 0.045);
                    Txt(0.16, 0.85, "black full: #sigma = 5.4 MeV; open: 6.2 MeV; teal: (a)", 0.033);
                    Txt(0.16, 0.80, "grey: sim unscaled, Gaussian at 0.5460, 6.2 MeV", 0.033, kGray + 2);
                    Txt(0.16, 0.75, "orange: #alpha_{peak} of template fit M1 (K0S scale)", 0.033, kOrange + 7);
                    continue;
                }
                struct V { int ir, sc, gv, mst, col; double dx; };
                const V vs[4] = {{0, 0, 0, 20, kBlack, 0.2}, {1, 0, 0, 24, kBlack, 0.4}, {0, 0, 1, 21, kColSim, 0.6}, {0, 1, 0, 22, kOrange + 7, 0.8}};
                for (const V& vv : vs)
                    for (int kk = 0; kk < (pans[p].k2 >= 0 ? 2 : 1); ++kk) {
                        TGraphErrors* g = new TGraphErrors();
                        for (int is = 0; is < kNSel; ++is) {
                            const Store& so = st[is][vv.ir][vv.sc][vv.gv][pans[p].im];
                            const int k = kk ? pans[p].k2 : pans[p].k;
                            const double v = k < 0 ? so.r.chi2 / so.r.ndf : so.r.a[k];
                            const double e = k < 0 ? 0 : so.r.ea[k];
                            g->SetPoint(is, is + vv.dx + (kk ? 0.07 : 0), v); g->SetPointError(is, 0, e);
                        }
                        int mst = vv.mst;
                        if (kk) mst = vv.mst == 20 ? 24 : vv.mst == 21 ? 25 : vv.mst == 22 ? 26 : 24;
                        if (pans[p].k2 >= 0 && vv.ir == 1) continue;   // keep two-parameter panels readable
                        g->SetMarkerStyle(mst); g->SetMarkerColor(vv.col); g->SetLineColor(vv.col); g->SetMarkerSize(1.4); g->SetLineWidth(2);
                        g->Draw("P same");
                    }
                Txt(0.15, 0.92, Form("(%c) %s", 'a' + p, pans[p].t), 0.045);
            }
            cv->SaveAs((fig + "/eta_enriched_fit_alpha.png").c_str());
            cv->SaveAs((fig + "/eta_enriched_fit_alpha.pdf").c_str());
        }
        printf("part 4: %.1f s\n", sw.RealTime());
    }

    // ======================================================================
    // part 8: subsystems in the peak and bump windows, E1 and E2
    // ======================================================================
    if (parts & 8) {
        TStopwatch sw;
        struct VarDef { const char* tag; const char* xt; const char* unit; int nb; double lo, hi; bool logy; };
        const int nV = 6;
        const VarDef vd[nV] = {{"mpipi", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "nb/(GeV/c^{2})", 20, 0.27, 0.47, false},
                               {"mee", "M_{e^{+}e^{-}} [GeV/c^{2}]", "nb/(GeV/c^{2})", 30, 0.0, 0.30, true},
                               {"oa", "OA_{e^{+}e^{-}} [deg]", "nb/deg", 30, 0.0, 60.0, true},
                               {"ppipi", "p_{#pi^{+}#pi^{-}} [GeV/c]", "nb/(GeV/c)", 25, 0.0, 2.5, false},
                               {"pee", "p_{e^{+}e^{-}} [GeV/c]", "nb/(GeV/c)", 25, 0.0, 1.5, false},
                               {"mm4", "MM_{4} [GeV/c^{2}]", "nb/(GeV/c^{2})", 30, 1.6, 3.1, false}};
        auto simVar = [](int iv) -> Var {
            switch (iv) {
                case 0: return [](const Row& r) { return (double)r.mpipi; };
                case 1: return [](const Row& r) { return (double)r.mee; };
                case 2: return [](const Row& r) { return (double)r.oa; };
                case 3: return [](const Row& r) { return r.p_pipi / 1000.; };
                case 4: return [](const Row& r) { return r.p_ee / 1000.; };
            }
            return [](const Row& r) { return (double)r.mm4; };
        };
        std::vector<DVar> dvars = {[](const scan::Row& r) { return (double)r.mpipi; }, [](const scan::Row& r) { return (double)r.mee; },
                                   [](const scan::Row& r) { return (double)r.oa; },
                                   [](const scan::Row& r) { return PairKinOf(r).p_pipi / 1000.; },
                                   [](const scan::Row& r) { return PairKinOf(r).p_ee / 1000.; }, [](const scan::Row& r) { return (double)r.mm4; }};
        std::vector<Binning> dbins;
        for (int iv = 0; iv < nV; ++iv) { Binning b; b.nbin = vd[iv].nb; b.lo = vd[iv].lo; b.hi = vd[iv].hi; dbins.push_back(b); }
        gTab = fopen((tab + "/eta_enriched_subsys.txt").c_str(), "w");
        P("# eta_enriched_subsys: subsystem distributions in the peak window (0.5336-0.5542) and the bump window (0.455-0.5225), selections E1 and E2\n");
        P("# data: N - CB per bin of the variable (CB = 2 sqrt(N++ N--) per bin; N with event-group variance); simulation T1 + T2 (and (a)), sum w / L_sim\n");
        P(Form("# absolute: data/sim content ratio in the histogram range (%s); shape: chi2 of the two normalised histograms (data + sim errors), ndf = bins - 1\n", kLdTxt));
        P("# note: the sum of per-bin CB differs from the CB of the window computed in M4 bins (the CB is not additive); both are listed\n");
        P(Form("%-4s %-5s %-6s %10s %10s %10s %10s %16s %16s %14s %14s %9s %9s %10s %10s\n", "sel", "win", "var", "N-CB rows", "CB rows", "CB(M4)", "outside",
               "data [nb]", "sim [nb]", "data/sim", "data/sim (a)", "chi2/ndf", "(a)", "mean data", "mean sim"));
        const int ssel[2] = {sE1, sE2};
        const char* wtag[2] = {"peak", "bump"};
        const double wlo[2] = {scan::kPeakLo, scan::kBumpLo}, whi[2] = {scan::kPeakHi, scan::kBumpHi};
        for (int js = 0; js < 2; ++js)
            for (int iw = 0; iw < 2; ++iw) {
                const int is = ssel[js];
                const std::vector<DH> dh = FillDH(dd, DatSelM(is, wlo[iw], whi[iw]), dvars, dbins);
                const scan::WinSum& wM4 = dReg[is].win[iw == 0 ? 2 + 3 : 2 + 1];
                std::vector<std::vector<WSpec>> sh(nV);
                for (int iv = 0; iv < nV; ++iv) {
                    Binning b = dbins[iv]; b.x = simVar(iv);
                    sh[iv] = FillMany(d, SimSelM(is, wlo[iw], whi[iw]), comps, opt, b);
                }
                TCanvas* cv = new TCanvas(Form("cSub_%d_%d", js, iw), "", 3300, 1600);
                cv->cd();
                MkPad(Form("pHsub%d%d", js, iw), 0, 0.86, 1, 1, false, 0, 0, 0, 0);
                Txt(0.005, 0.80, Form("Subsystems in the %s window (M_{4} %.4f-%.4f GeV/c^{2}), %s: data N - CB against the simulation; top absolute, bottom normalised to the content",
                                      iw == 0 ? "peak" : "bump", wlo[iw], whi[iw], SelLatex(is).c_str()), 0.13);
                Txt(0.005, 0.60, NormLine().c_str(), 0.085);
                Txt(0.005, 0.46, Form("%s; data CB = 2#sqrt{N_{++}N_{--}} per bin of the variable; shape #chi^{2}/ndf: normalised data vs normalised T1 + T2", PreLine()), 0.085);
                SharedLegend(0.0, 0.0, 1.0, 0.40, 6, 0.08, true, false, false);
                for (int iv = 0; iv < nV; ++iv) {
                    const DH& h = dh[iv];
                    const double bw = h.bw;
                    TH1D* hD = h.Hist(Form("sD_%d%d%d", js, iw, iv), 2, 1. / (kLdata * bw));
                    TH1D* hSum = sh[iv][cSum].Hist(Form("sS_%d%d%d", js, iw, iv), 2);
                    TH1D* hSumA = sh[iv][cSumA].Hist(Form("sA_%d%d%d", js, iw, iv), 2);
                    const double Dn = h.Sum(true), Dv = h.SumVar();
                    const double Sn = sh[iv][cSum].S, Sv = sh[iv][cSum].V, An = sh[iv][cSumA].S;
                    const double rr = (Dn / kLdata) / (Sn / kLsim);
                    const double rre = rr * std::sqrt(Dv / (Dn * Dn) + Sv / (Sn * Sn));
                    const double ra = (Dn / kLdata) / (An / kLsim);
                    // shape chi2
                    auto shapeChi2 = [&](const WSpec& s, double& chi, int& ndf) {
                        chi = 0; ndf = -1;
                        for (int i = 0; i < h.nb; ++i) {
                            const double a = h.sig(i) / Dn, ea2 = h.sigVar(i) / (Dn * Dn);
                            const double bq = s.s[i] / s.S, eb2 = s.v[i] / (s.S * s.S);
                            if (ea2 + eb2 <= 0 || (h.n[i] == 0 && s.s[i] == 0)) continue;
                            chi += (a - bq) * (a - bq) / (ea2 + eb2); ++ndf;
                        }
                    };
                    double chi, chiA; int ndf, ndfA;
                    shapeChi2(sh[iv][cSum], chi, ndf); shapeChi2(sh[iv][cSumA], chiA, ndfA);
                    double md = 0, ms = 0;
                    for (int i = 0; i < h.nb; ++i) { const double x = h.lo + (i + 0.5) * bw; md += x * h.sig(i); ms += x * sh[iv][cSum].s[i]; }
                    md /= Dn; ms /= Sn;
                    double cbs = 0; for (int i = 0; i < h.nb; ++i) cbs += h.cb(i);
                    // rows outside the range (data N, simulation weight share)
                    const double outS = sh[iv][cSum].Sall > 0 ? 1. - sh[iv][cSum].S / sh[iv][cSum].Sall : 0;
                    P(Form("%-4s %-5s %-6s %10.0f %10.0f %10.0f %9.4f%% %7.4f+-%-7.4f %7.4f+-%-7.4f %6.3f+-%-6.3f %6.3f%8s %6.1f/%-3d %6.1f/%-3d %10.4f %10.4f\n",
                           kSel[is].tag, wtag[iw], vd[iv].tag, Dn, cbs, wM4.cb, 100 * outS, Dn / kLdata, std::sqrt(Dv) / kLdata, Sn / kLsim,
                           std::sqrt(Sv) / kLsim, rr, rre, ra, "", chi, ndf, chiA, ndfA, md, ms));
                    for (int row = 0; row < 2; ++row) {
                        const bool norm = row == 1;
                        cv->cd();
                        const double x1 = iv / 6., x2 = (iv + 1) / 6., y1 = norm ? 0.0 : 0.43, y2 = norm ? 0.43 : 0.86;
                        MkPad(Form("ps_%d%d%d%d", js, iw, iv, row), x1, y1, x2, y2, vd[iv].logy, 0.17, 0.03, 0.12, 0.14);
                        const double fD = norm ? kLdata * bw / Dn : 1., fS = norm ? kLsim * bw / Sn : 1.;
                        TH1D* d1 = (TH1D*)hD->Clone(Form("d1_%d%d%d%d", js, iw, iv, row)); d1->Scale(fD);
                        TH1D* s1 = (TH1D*)hSum->Clone(Form("s1_%d%d%d%d", js, iw, iv, row)); s1->Scale(fS);
                        TH1D* a1 = (TH1D*)hSumA->Clone(Form("a1_%d%d%d%d", js, iw, iv, row)); a1->Scale(norm ? kLsim * bw / An : 1.);
                        TH1D* hs[6];
                        for (int k = 0; k < 6; ++k) {
                            hs[k] = sh[iv][kStack[k]].Hist(Form("ss_%d%d%d%d%d", js, iw, iv, row, k), 2);
                            hs[k]->Scale(fS);
                            if (k > 0) hs[k]->Add(hs[k - 1]);
                            StyleH(hs[k], kStackLine[k], 1, 1, kStackFill[k]);
                        }
                        double ym = std::max(HMax(d1), std::max(HMax(s1), HMax(a1)));
                        double ymin = std::max(1e-6, std::min(HMinPos(s1), HMinPos(d1)) * 0.3);
                        TH1D* fr = (TH1D*)d1->Clone(Form("sfr_%d%d%d%d", js, iw, iv, row)); fr->Reset();
                        Axes(fr, vd[iv].xt, norm ? "fraction per bin" : Form("d#sigma/dx [%s]", vd[iv].unit), 0.060, 1.35);
                        if (vd[iv].logy) { fr->SetMinimum(std::max(ymin, ym * 1e-4)); fr->SetMaximum(ym * 4); }
                        else { fr->SetMinimum(0); fr->SetMaximum(ym * 1.25); }
                        fr->Draw("axis");
                        if (vd[iv].logy) {
                            for (int k = 5; k >= 0; --k) { hs[k]->SetFillStyle(0); hs[k]->SetLineWidth(2); }
                            TH1D* hk = sh[iv][cKM].Hist(Form("skm_%d%d%d%d", js, iw, iv, row), 2); hk->Scale(fS);
                            StyleH(hk, kColPhysBkg, 1, 1, kColPhysFill); hk->Draw("hist same");
                            for (int k = 1; k < 6; ++k) {
                                TH1D* hc = sh[iv][kStack[k]].Hist(Form("sc_%d%d%d%d%d", js, iw, iv, row, k), 2); hc->Scale(fS);
                                StyleH(hc, kStackLine[k], 3, 1); hc->Draw("hist same");
                            }
                        } else
                            for (int k = 5; k >= 0; --k) hs[k]->Draw("hist same");
                        StyleH(s1, kColSim, 4, 1); s1->Draw("hist same");
                        StyleH(a1, kColSim, 3, 2); a1->Draw("hist same");
                        StyleM(d1, kColSig, 20, 1.1); d1->Draw("E1 same");
                        fr->Draw("axis same");
                        if (!norm)
                            Txt(0.19, 0.91, Form("(%c) data/sim %.2f #pm %.2f, (a) %.2f", 'a' + iv, rr, rre, ra), 0.050);
                        else
                            Txt(0.19, 0.91, Form("(%c') shape #chi^{2}/ndf %.0f/%d, (a) %.0f/%d", 'a' + iv, chi, ndf, chiA, ndfA), 0.048);
                        if (!vd[iv].logy && !norm && iv == 0)
                            Txt(0.19, 0.84, "linear panels: stacked; log panels: lines", 0.045);
                    }
                }
                const std::string fn = fig + Form("/eta_enriched_sub_%s_%s", kSel[is].tag, wtag[iw]);
                cv->SaveAs((fn + ".png").c_str());
                if (js == 0 || iw == 0) cv->SaveAs((fn + ".pdf").c_str());
            }
        fclose(gTab); gTab = nullptr;
        printf("part 8: %.1f s\n", sw.RealTime());
    }

    fprintf(fchk, "MM4 cut used for E1mm, E2mm: %.2f GeV/c2\n", gMMcut);
    fclose(fchk);
    printf("sim_eta_enriched: total %.1f s\n", swAll.RealTime());
}
