// ============================================================================
// wp1_repro.C  --  WP1 "before": reproduction of the simulation-study table
//                  rw_compare_enriched.txt (selection pre) and the nominal pre/S0
//                  simulated window sums at s = 1 and s = 0.99615 (prefix wp1)
// ============================================================================
// Why: hists_v1.root holds the simulation only with the mass transfer of D8
// (M4' = 0.99615 M4). The simulation study compared the data with the UNSCALED
// simulated M4 (s = 1), L_data = 6080.5 nb^-1 and no 64/63. To test the
// acceptance criterion of WP1 (pre reproduces rw_compare_enriched.txt x 1.0161
// within 0.5 %) the simulated window sums are recomputed here at s = 1 AND at
// s = 0.99615, by streaming the simulated rows (no full load, < 1 GB):
//   selections pre and S0 (defs.h), mode A, w x physics code 0 (before) and 108
//   (after, V3), T1 + T2 and the ghost variant (a), event-group errors.
// The s = 0.99615 sums must equal hists_v1.root (closure check, printed).
// The data are read from hists_v1.root (pre, raw rows; windows of the 'win'
// layout and the fine 2.5 MeV M4 spectrum).
// Output: docs/latex/tables_all/wp1_repro.txt
// Run (from /home/przygoda/HADES/FAT/FAT; about 1 min, < 1 GB):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp1", kTRUE)' 'plots/all_report/wp1_repro.C+O'
// ============================================================================
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "defs.h"
#include "hists_io.h"
#include "sim_eta_physics.h"

namespace wp1r {
using namespace alldefs;
const char* const kTab = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all";
const char* const kOldFile = "/home/przygoda/HADES/FAT/FAT_simwt/plots/eta_sim/tables/rw_compare_enriched.txt";
const double kLold = 6080.5;   // L_data of the simulation study (no 64/63)

// windows (data scale; simulated rows use M4' = s M4)
struct W { const char* tag; double lo, hi; };
const int kNW = 20;
const W kW[kNW] = {
    {"below bump 0.35-0.455", 0.35, 0.455},          // 0  old regions (sim_eta_enriched.C kReg)
    {"bump window 0.455-0.5225", 0.4550, 0.5225},     // 1
    {"gap 0.5225-0.5336", 0.5225, 0.53364},           // 2
    {"peak window 0.5336-0.5542", 0.53364, 0.55416},  // 3
    {"0.5542-0.56", 0.55416, 0.56},                   // 4
    {"0.56-0.70", 0.56, 0.70},                        // 5
    {"0.70-1.00", 0.70, 1.00},                        // 6
    {"all 0.35-1.00", 0.35, 1.00},                    // 7
    {"band 0.4875-0.5225", 0.4875, 0.5225},           // 8  jump bands
    {"band 0.5600-0.5950", 0.56, 0.595},              // 9
    {"R0 <0.455", -1e9, 0.455},                       // 10 Table B1 regions
    {"R1 B", 0.4550, 0.5225},
    {"R2 0.5225-0.5336", 0.5225, 0.53364},
    {"R3 P", 0.53364, 0.55416},
    {"R4 0.5542-0.565", 0.55416, 0.565},
    {"R5 0.565-0.70", 0.565, 0.70},
    {"R6 0.70-0.90", 0.70, 0.90},
    {"R7 0.90-1.10", 0.90, 1.10},                     // 17
    {"Fit 0.40-0.70", 0.40, 0.70},                    // 18
    {"Bef 0.30-1.10", 0.30, 1.10}};                   // 19
// combinations: sel (0 pre, 1 S0) x phys (0 off, 1 = 108) x var (0 T12, 1 ghost (a)) x scale (0: s = 1, 1: s = 0.99615)
inline int Combo(int sel, int ph, int v, int sc) { return ((sel * 2 + ph) * 2 + v) * 2 + sc; }
const int kNC = 16;
inline int Q(int c, int w) { return c * kNW + w; }
const int kNQ = kNC * kNW;

struct Acc {
    std::vector<double> sum, var, grp;
    std::vector<int> touched;
    std::vector<char> isT;
    Acc() : sum(kNQ, 0.), var(kNQ, 0.), grp(kNQ, 0.), isT(kNQ, 0) {}
    void Add(int q, double w) {
        if (!isT[q]) { isT[q] = 1; touched.push_back(q); }
        grp[q] += w;
        sum[q] += w;
    }
    void Flush() {
        for (int q : touched) { var[q] += grp[q] * grp[q]; grp[q] = 0.; isT[q] = 0; }
        touched.clear();
    }
};

struct VE { double v = 0, e = 0; };
VE Ratio(VE a, VE b) {
    VE r;
    if (!(b.v != 0)) return r;
    r.v = a.v / b.v;
    r.e = std::fabs(r.v) * std::sqrt(std::pow(a.e / a.v, 2) + std::pow(b.e / b.v, 2));
    return r;
}
double RepRms(const TH1* rep, int bin, int nrep) {
    double m = 0, m2 = 0;
    for (int r = 1; r <= nrep; ++r) { const double x = rep->GetBinContent(bin, r); m += x; m2 += x * x; }
    m /= nrep;
    return std::sqrt(std::max(0., m2 / nrep - m * m));
}

// old-table values for "pre": data N-CB [nb], sim before/after [nb], ratio before (a), after (a)
struct Old { bool ok = false; double d = 0, de = 0, sb = 0, sa = 0, rb = 0, rbe = 0, rba = 0, ra = 0, rae = 0, raa = 0; };
std::map<std::string, Old> ReadOld() {
    std::map<std::string, Old> m;
    std::ifstream in(kOldFile);
    std::string line;
    bool sec1 = false;
    while (std::getline(in, line)) {
        if (line.rfind("## 1.", 0) == 0) { sec1 = true; continue; }
        if (line.rfind("## 2.", 0) == 0) break;
        if (!sec1 || line.rfind("pre ", 0) != 0) continue;
        const char* names[] = {"below bump 0.35-0.455", "bump window 0.455-0.5225", "gap 0.5225-0.5336", "peak window 0.5336-0.5542", "0.5542-0.56",
                               "0.56-0.70", "0.70-1.00", "all 0.35-1.00", "band 0.4875-0.5225", "band 0.5600-0.5950", "JUMP = band 0.56 / band 0.49"};
        for (const char* nm : names) {
            const size_t p = line.find(nm);
            if (p == std::string::npos) continue;
            Old o;
            const size_t bar = line.find('|');
            if (std::string(nm).rfind("band", 0) == 0 || std::string(nm).rfind("JUMP", 0) == 0) {
                o.ok = sscanf(line.c_str() + bar + 1, " %lf +- %lf (%lf) %lf +- %lf (%lf)", &o.rb, &o.rbe, &o.rba, &o.ra, &o.rae, &o.raa) == 6;
            } else {
                const char* s = line.c_str() + p + strlen(nm);
                const int n1 = sscanf(s, " %lf+-%lf %lf -> %lf", &o.d, &o.de, &o.sb, &o.sa);
                const int n2 = sscanf(line.c_str() + bar + 1, " %lf +- %lf (%lf) %lf +- %lf (%lf)", &o.rb, &o.rbe, &o.rba, &o.ra, &o.rae, &o.raa);
                o.ok = (n1 == 4 && n2 == 6);
            }
            m[nm] = o;
            break;
        }
    }
    return m;
}
}  // namespace wp1r

void wp1_repro() {
    using namespace wp1r;
    TStopwatch sw;
    // ---------------------------------------------------------------- simulation, streaming
    simphys::Init();
    simphys::LoadAngles();
    Acc A;
    struct TF { const char* file; const char* tree; int tid; };
    const TF tf[3] = {{simeng::kCompact, "std_us", 0}, {simeng::kCompact, "lep_us", 3}, {simeng::kCompactEta, "eta_us", 6}};
    long long nGrpBack = 0, nUsed = 0;
    for (const TF& x : tf) {
        TFile* f = TFile::Open(x.file, "READ");
        TTree* t = (TTree*)f->Get(x.tree);
        simeng::Row r;
        simeng::BindRow(t, r, false);
        const Long64_t n = t->GetEntries();
        std::unordered_set<long long> done;
        long long cur = -1;
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            const long long gk = ((long long)x.tid << 40) | (long long)(unsigned)r.grp;
            if (gk != cur) {
                A.Flush();
                if (cur >= 0) done.insert(cur);
                if (done.count(gk)) ++nGrpBack;
                cur = gk;
            }
            if (!r.ownA) continue;                       // mode A
            if (!simeng::PassPre(r)) continue;
            const bool s0 = PiVtx(r);
            const bool x17 = (r.tmpl == 1 && r.cls == 6);
            if (x17) continue;                           // X17 never in a sum
            const bool t12 = (r.tmpl == 1 && r.cls >= 2 && r.cls <= 5) || r.tmpl == 2;
            const bool ga = t12 || simeng::GhostPair(r);
            if (!ga) continue;
            ++nUsed;
            const double w0 = r.w, w1 = r.w * simphys::RowFactor(r, kPhysNom);
            for (int sc = 0; sc < 2; ++sc) {
                const double m = sc == 0 ? (double)r.m4 : M4sim(r.m4, r.m4_true, kS, kR);
                for (int wi = 0; wi < kNW; ++wi) {
                    if (!(m >= kW[wi].lo && m < kW[wi].hi)) continue;
                    for (int sel = 0; sel < 2; ++sel) {
                        if (sel == 1 && !s0) continue;
                        for (int ph = 0; ph < 2; ++ph) {
                            const double w = ph ? w1 : w0;
                            if (t12) A.Add(Q(Combo(sel, ph, 0, sc), wi), w);
                            A.Add(Q(Combo(sel, ph, 1, sc), wi), w);
                        }
                    }
                }
            }
        }
        A.Flush();
        printf("wp1_repro: %s %lld rows\n", x.tree, (long long)n);
        t->ResetBranchAddresses();
        f->Close();
    }
    printf("wp1_repro: %lld selected owned rows used; event groups that reappear non-adjacently: %lld\n", nUsed, nGrpBack);
    auto Sim = [&](int sel, int ph, int v, int sc, int wi) {
        VE r; const int q = Q(Combo(sel, ph, v, sc), wi);
        r.v = A.sum[q] / kLsim; r.e = std::sqrt(A.var[q]) / kLsim; return r;
    };

    // ---------------------------------------------------------------- data from hists_v1.root
    TFile* hf = TFile::Open(kHistFile);
    FILE* o = fopen(Form("%s/wp1_repro.txt", kTab), "w");
    fprintf(o, "# wp1_repro.txt -- plots/all_report/wp1_repro.C (WP1 of SPEC report_all)\n");
    fprintf(o, "# Reproduction of FAT_simwt/plots/eta_sim/tables/rw_compare_enriched.txt (selection pre, section 1) and the effect of the\n");
    fprintf(o, "# mass transfer D8. Simulation streamed here (std + lep + eta owned rows, mode A, w x code 0 'before' / 108 'after', T1 + T2;\n");
    fprintf(o, "# (a) = + ghost-flagged true pairs), sigma = sum w / %.2f nb^-1, event-group errors. Data N - CB from hists_v1.root (pre, raw rows,\n", kLsim);
    fprintf(o, "# CB = spread rule; the B, P and jump-band windows are exact 'win' cells with CB replica errors; the other old regions are sums of\n");
    fprintf(o, "# the fine 2.5 MeV bins (0.5542-0.56 = R4 minus the fine bins 0.560-0.565), errors in quadrature).\n");
    fprintf(o, "# selected owned rows used: %lld; non-adjacent event groups: %lld\n\n", nUsed, nGrpBack);

    TH1* nw = allio::Data(hf, "pre", "all", "N", "win");
    TH1* cw = allio::Data(hf, "pre", "all", "CB", "win");
    TH2* cr = (TH2*)allio::Data(hf, "pre", "all", "CBrep", "win");
    TH1* nf = allio::Data(hf, "pre", "all", "N", "fine_M4");
    TH1* cf = allio::Data(hf, "pre", "all", "CB", "fine_M4");
    auto Fine = [&](double a, double b, double& n, double& nv, double& c, double& cv) {
        n = nv = c = cv = 0;
        for (int i = 1; i <= nf->GetNbinsX(); ++i) {
            const double x = nf->GetXaxis()->GetBinCenter(i);
            if (x < a || x > b) continue;
            n += nf->GetBinContent(i); nv += std::pow(nf->GetBinError(i), 2);
            c += cf->GetBinContent(i); cv += std::pow(cf->GetBinError(i), 2);
        }
    };
    // data N - CB in raw rows for the old regions 0..9
    VE draw[10];
    for (int k = 0; k < 10; ++k) {
        int wid = -1;
        if (k == 1) wid = kWB; else if (k == 3) wid = kWP; else if (k == 8) wid = kWJlo; else if (k == 9) wid = kWJhi; else if (k == 2) wid = kWR2;
        if (wid >= 0) {
            draw[k].v = allio::Win(nw, wid) - allio::Win(cw, wid);
            draw[k].e = std::hypot(allio::WinErr(nw, wid), RepRms(cr, wid + 1, kNRepCB));
            continue;
        }
        double n, nv, c, cv;
        if (k == 4) {
            double n2, nv2, c2, cv2;
            Fine(0.56, 0.565, n2, nv2, c2, cv2);
            n = allio::Win(nw, kWR4) - n2; c = allio::Win(cw, kWR4) - c2;
            nv = std::max(0., std::pow(allio::WinErr(nw, kWR4), 2) - nv2); cv = std::pow(RepRms(cr, kWR4 + 1, kNRepCB), 2);
        } else {
            Fine(kW[k].lo, kW[k].hi, n, nv, c, cv);
        }
        draw[k].v = n - c; draw[k].e = std::sqrt(nv + cv);
    }
    // jump bands use the 'win' cells; J = R(hi band) / R(lo band)
    std::map<std::string, Old> old = ReadOld();
    const double fOld = kPT3 * kLold / kL;   // 1.01609: factor between this normalisation and the old one

    fprintf(o, "## 1. Reproduction at the settings of the simulation study (s = 1, L_data = 6080.5 nb^-1, no 64/63), pre\n");
    fprintf(o, "%-27s | %-17s %-17s | %-17s %-17s | %-29s %-29s | %s\n", "region", "data N-CB [nb]", "old", "sim bef -> aft", "old", "data/sim before (a)",
            "data/sim after (a)", "old after (a) ; rel. diff of 'after'");
    double maxDiff = 0;
    VE rs1[10][2][2];   // [region][phys][var] at s = 1
    for (int k = 0; k < 10; ++k) {
        const VE d = {draw[k].v / kLold, draw[k].e / kLold};
        for (int ph = 0; ph < 2; ++ph)
            for (int v = 0; v < 2; ++v) rs1[k][ph][v] = Ratio(d, Sim(0, ph, v, 0, k));
        const Old& x = old[kW[k].tag];
        const double rel = x.ok ? rs1[k][1][0].v / x.ra - 1. : NAN;
        if (x.ok && std::isfinite(rel)) maxDiff = std::max(maxDiff, std::fabs(rel));
        const std::string od = k >= 8 ? std::string("   -    ") : std::string(Form("%8.4f", x.d));
        const std::string os = k >= 8 ? std::string("        -        ") : std::string(Form("%7.3f -> %-7.3f", x.sb, x.sa));
        fprintf(o, "%-27s | %8.4f+-%-7.4f %s          | %7.3f -> %-7.3f %s | %6.3f +- %-5.3f (%5.3f)      %6.3f +- %-5.3f (%5.3f)      | %6.3f (%5.3f) ; %+.2e\n",
                kW[k].tag, d.v, d.e, od.c_str(), Sim(0, 0, 0, 0, k).v, Sim(0, 1, 0, 0, k).v, os.c_str(), rs1[k][0][0].v, rs1[k][0][0].e, rs1[k][0][1].v,
                rs1[k][1][0].v, rs1[k][1][0].e, rs1[k][1][1].v, x.ra, x.raa, rel);
    }
    VE J1[2][2];
    for (int ph = 0; ph < 2; ++ph)
        for (int v = 0; v < 2; ++v) J1[ph][v] = Ratio(rs1[9][ph][v], rs1[8][ph][v]);
    {
        const Old& x = old["JUMP = band 0.56 / band 0.49"];
        fprintf(o, "%-27s | %35s | %35s | %6.3f +- %-5.3f (%5.3f)      %6.3f +- %-5.3f (%5.3f)      | %6.3f (%5.3f) ; %+.2e\n", "JUMP = band 0.56 / band 0.49", "", "",
                J1[0][0].v, J1[0][0].e, J1[0][1].v, J1[1][0].v, J1[1][0].e, J1[1][1].v, x.ra, x.raa, J1[1][0].v / x.ra - 1.);
        maxDiff = std::max(maxDiff, std::fabs(J1[1][0].v / x.ra - 1.));
    }
    fprintf(o, "# largest relative difference of the 'after' ratio to the old table (regions + jump): %.2e\n\n", maxDiff);

    fprintf(o, "## 2. Nominal normalisation of this report (L = %.1f nb^-1, x 64/63) and the mass transfer D8, pre, after (code 108)\n", kL);
    fprintf(o, "# old x 1.0161 = the old 'after' ratio with this normalisation; s = 1 / s = 0.99615 = the ratio with the simulated M4 unscaled / scaled\n");
    fprintf(o, "%-27s | %-12s | %-17s | %-17s | %-17s | %-9s | %-9s\n", "region", "old x 1.0161", "s = 1", "s = 0.99615 (nom)", "(a), s = 0.99615",
            "s=1 / old", "nom / old");
    VE rn[10][2];
    for (int k = 0; k < 10; ++k) {
        const VE d = {draw[k].v * kRawToNb, draw[k].e * kRawToNb};
        const VE r1 = Ratio(d, Sim(0, 1, 0, 0, k));
        rn[k][0] = Ratio(d, Sim(0, 1, 0, 1, k));
        rn[k][1] = Ratio(d, Sim(0, 1, 1, 1, k));
        const Old& x = old[kW[k].tag];
        fprintf(o, "%-27s | %12.4f | %7.4f +- %-6.4f | %7.4f +- %-6.4f | %7.4f +- %-6.4f | %+8.2f%% | %+8.2f%%\n", kW[k].tag, x.ra * kOldTab, r1.v, r1.e, rn[k][0].v,
                rn[k][0].e, rn[k][1].v, rn[k][1].e, 100. * (r1.v / (x.ra * kOldTab) - 1.), 100. * (rn[k][0].v / (x.ra * kOldTab) - 1.));
    }
    const VE Jn = Ratio(rn[9][0], rn[8][0]), Jna = Ratio(rn[9][1], rn[8][1]);
    {
        const Old& x = old["JUMP = band 0.56 / band 0.49"];
        const VE d8 = {draw[8].v * kRawToNb, draw[8].e * kRawToNb}, d9 = {draw[9].v * kRawToNb, draw[9].e * kRawToNb};
        const VE J1n = Ratio(Ratio(d9, Sim(0, 1, 0, 0, 9)), Ratio(d8, Sim(0, 1, 0, 0, 8)));
        fprintf(o, "%-27s | %12.4f | %7.4f +- %-6.4f | %7.4f +- %-6.4f | %7.4f +- %-6.4f | %+8.2f%% | %+8.2f%%\n", "JUMP (no normalisation)", x.ra, J1n.v, J1n.e, Jn.v, Jn.e,
                Jna.v, Jna.e, 100. * (J1n.v / x.ra - 1.), 100. * (Jn.v / x.ra - 1.));
    }
    fprintf(o, "# factor of the normalisation (64/63)(6080.5/6079.2) = %.5f (defs.h kOldTab %.4f)\n\n", fOld, kOldTab);

    // ---------------------------------------------------------------- closure against hists_v1.root (s = 0.99615)
    fprintf(o, "## 3. Closure: streamed s = 0.99615 sums against hists_v1.root (sim/<sel>__<var>/T12_win), Table B1 regions R0-R7, Fit, Bef\n");
    const char* selT[2] = {"pre", "S0"};
    const char* varT[3] = {"wtsOff", "nom", "ghostA"};
    int nBad = 0, nChk = 0;
    double worst = 0;
    for (int sel = 0; sel < 2; ++sel)
        for (int vv = 0; vv < 3; ++vv) {
            const int ph = vv == 0 ? 0 : 1, v = vv == 2 ? 1 : 0;
            TH1* h = allio::Sim(hf, selT[sel], varT[vv], kT12, "win");
            if (!h) continue;
            const int wid[10] = {kWR0, kWR1, kWR2, kWR3, kWR4, kWR5, kWR6, kWR7, kWFit, kWBef};
            for (int j = 0; j < 10; ++j) {
                const VE s = Sim(sel, ph, v, 1, 10 + j);
                const double ref = allio::Win(h, wid[j]), rel = ref != 0 ? s.v / ref - 1. : 0.;
                const double rele = allio::WinErr(h, wid[j]) > 0 ? s.e / allio::WinErr(h, wid[j]) - 1. : 0.;
                ++nChk;
                worst = std::max(worst, std::fabs(rel));
                if (std::fabs(rel) > 1e-9 || std::fabs(rele) > 1e-6) {
                    ++nBad;
                    fprintf(o, "MISMATCH %s__%s %s: streamed %.6f +- %.6f, file %.6f +- %.6f\n", selT[sel], varT[vv], kW[10 + j].tag, s.v, s.e, ref,
                            allio::WinErr(h, wid[j]));
                }
            }
        }
    fprintf(o, "# %d comparisons, %d mismatches, largest relative difference of the content %.2e\n\n", nChk, nBad, worst);

    // ---------------------------------------------------------------- acceptance
    const Old& ob = old["bump window 0.455-0.5225"];
    const Old& op = old["peak window 0.5336-0.5542"];
    const Old& oj = old["JUMP = band 0.56 / band 0.49"];
    fprintf(o, "## 4. Acceptance criterion of WP1: pre reproduces rw_compare_enriched.txt x 1.0161 within 0.5 %% (bump 0.958, peak 1.233); J 1.44-1.47\n");
    auto pf = [](bool b) { return b ? "PASS" : "FAIL"; };
    const double b1 = Ratio({draw[1].v * kRawToNb, 0}, Sim(0, 1, 0, 0, 1)).v, p1 = Ratio({draw[3].v * kRawToNb, 0}, Sim(0, 1, 0, 0, 3)).v;
    fprintf(o, "(i)  same simulation settings as the old table (s = 1), this normalisation: bump %.4f vs %.4f (%+.2f%%) %s; peak %.4f vs %.4f (%+.2f%%) %s\n", b1,
            ob.ra * kOldTab, 100 * (b1 / (ob.ra * kOldTab) - 1), pf(std::fabs(b1 / (ob.ra * kOldTab) - 1) < 0.005), p1, op.ra * kOldTab,
            100 * (p1 / (op.ra * kOldTab) - 1), pf(std::fabs(p1 / (op.ra * kOldTab) - 1) < 0.005));
    fprintf(o, "     J (s = 1) = %.4f +- %.4f %s (range 1.44-1.47)\n", J1[1][0].v, J1[1][0].e, pf(J1[1][0].v >= 1.44 && J1[1][0].v <= 1.47));
    fprintf(o, "(ii) nominal (s = 0.99615, D8): bump %.4f +- %.4f vs %.4f (%+.2f%%) %s; peak %.4f +- %.4f vs %.4f (%+.2f%%) %s\n", rn[1][0].v, rn[1][0].e,
            ob.ra * kOldTab, 100 * (rn[1][0].v / (ob.ra * kOldTab) - 1), pf(std::fabs(rn[1][0].v / (ob.ra * kOldTab) - 1) < 0.005), rn[3][0].v, rn[3][0].e,
            op.ra * kOldTab, 100 * (rn[3][0].v / (op.ra * kOldTab) - 1), pf(std::fabs(rn[3][0].v / (op.ra * kOldTab) - 1) < 0.005));
    fprintf(o, "     J (nominal) = %.4f +- %.4f %s (range 1.44-1.47); old %.3f\n", Jn.v, Jn.e, pf(Jn.v >= 1.44 && Jn.v <= 1.47), oj.ra);
    fprintf(o, "# (i) tests the chain (data, CB, simulation, weights, units); (ii) adds the deliberate mass transfer of D8, which moves simulated rows\n");
    fprintf(o, "#     across the window edges (sim content s = 1 -> 0.99615: bump %.3f -> %.3f nb, peak %.3f -> %.3f nb, jump bands %.3f -> %.3f / %.3f -> %.3f nb)\n",
            Sim(0, 1, 0, 0, 1).v, Sim(0, 1, 0, 1, 1).v, Sim(0, 1, 0, 0, 3).v, Sim(0, 1, 0, 1, 3).v, Sim(0, 1, 0, 0, 8).v, Sim(0, 1, 0, 1, 8).v, Sim(0, 1, 0, 0, 9).v,
            Sim(0, 1, 0, 1, 9).v);

    // ---------------------------------------------------------------- S0 at both scales (for the results note)
    fprintf(o, "\n## 5. S0 (nominal selection), after, T1 + T2: simulated nb at s = 1 and s = 0.99615 in the old regions\n");
    for (int k = 0; k < 10; ++k)
        fprintf(o, "%-27s | s = 1 %9.4f +- %-7.4f | s = 0.99615 %9.4f +- %-7.4f | (a) %9.4f\n", kW[k].tag, Sim(1, 1, 0, 0, k).v, Sim(1, 1, 0, 0, k).e, Sim(1, 1, 0, 1, k).v,
                Sim(1, 1, 0, 1, k).e, Sim(1, 1, 1, 1, k).v);
    fclose(o);
    printf("wp1_repro: wrote %s/wp1_repro.txt in %.1f s\n", kTab, sw.RealTime());
}
