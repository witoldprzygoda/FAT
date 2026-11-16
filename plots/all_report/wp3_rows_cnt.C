// wp3_rows_cnt.C -- COPY of wp3_rows.C (figure agent B, 2026-11-15): "counts" instead of "rows" in the y title and the
// step-height annotation of Fig. wp3_crk_edge (its two annotation lines start at x = 0.175 instead of 0.19 so that the longer
// word keeps the line inside the frame); writes ONLY docs/latex/figs_all/wp3_crk_edge_cnt.{pdf,png}; wp3_rows.root/.txt,
// the summary and wp3_crk_tag go to /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/wp3_rows/ (the WP3 files are not touched).
// Run: root -l -b -q -e 'gSystem->SetBuildDir("<scratch>/aclic_wp3_cnt", kTRUE)' 'plots/all_report/wp3_rows_cnt.C+O'
// ============================================================================
// wp3_rows.C  --  WP3 (control regions, stage A): the parts that need the data ROWS
// ============================================================================
// SPEC report_all, WP3. Written by WP3 (prefix wp3), 2026-11-14.
//
//   1. Event-group covariance of the data N over the 700 F10 cells (S0 and pre):
//      Sigma_ab = sum over event groups g of n_ga n_gb (n_ga = selected rows of group g in cell a).
//      Rows are grouped by ADJACENCY of grp (never by a global map). Used by wp3_stageA.C for the
//      sandwich errors. Check: the per-cell counts equal N_F10 of hists_v1.root exactly.
//   2. Test A2, CR-K null test (edge step). N - CB in the M(pi+pi-) slices 0.43-0.46 and 0.46-0.55
//      (and 0.44-0.46 for the fall-back "CR-K at 0.44"), fitted with pol2 + h x step at the FIXED edge
//      E(0.4116, 0) = 0.4739 GeV/c2 (data scale), smearing 5.4 MeV, range edge -+ 45 MeV (engine
//      scan::FitEdgeStep, the procedure of E_step). Two CB constructions:
//        "own"    = 2 sqrt(N++ N--) of the slice itself in 10 MeV cells (as the old E_step control),
//        "spread" = the report's rule: CB of all slices x f(slice | 50 MeV band) from N++ + N--.
//      Validation: the old control slice 0.42-0.46 on "none" with the hand-set edge 0.5150 must
//      give 36.5 +- 19.3 (scan_bump.C).
//   3. Test A2, second part: ECAL-tagged eta (frozen definition of scan_bump.C) in CR-K: m5c peak
//      content T with the frozen shape 0.5463 / 6.2 MeV, pol2, 0.49-0.62 (P_main). Positive control:
//      the same fit for tagged rows of the bump region (S1-S4, M4 < 0.5225).
//
// Outputs: /hdd3/przygoda/HADES/FAT_out/all_report/wp3/wp3_rows.root, wp3_rows.txt;
//          docs/latex/figs_all/wp3_crk_edge.pdf (+png), wp3_crk_tag.pdf (+png).
// Run (from /home/przygoda/HADES/FAT/FAT; 2.6 GB, about 1 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp3", kTRUE)' \
//        'plots/all_report/wp3_rows.C+O'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMatrixDSym.h"
#include "TStyle.h"
#include "TSystem.h"

#include "hists_io.h"
#include "s5_blind_ana_fit.h"

namespace wp3r {
using namespace alldefs;

const char* const kOutDir = "/tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figB/out/wp3_rows";   // _cnt copy: all non-figure outputs to scratch
const char* const kFigDir = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";

FILE* gOut = nullptr;
FILE* gSum = nullptr;   // machine-readable summary (read by wp3_stageA.C)
#define P(...) do { printf(__VA_ARGS__); if (gOut) fprintf(gOut, __VA_ARGS__); } while (0)

// ---- ECAL tag, copied from plots/scan_bump.C (frozen; bump_pi0.txt) ----------------------------
const double kPi0Lo = 0.1190, kPi0Hi = 0.1730;   // pi0 window in M(e+e-gamma)
const double kMeeMax = 0.100;
const double kTagMean = 0.5463, kTagSigma = 0.0062, kTagLo = 0.49, kTagHi = 0.62;
const int kTagOrder = 2, kTagNbin = 160, kTagReb = 8;
const double kTagHLo = 0.40, kTagHHi = 0.80;
inline double M5c(const scan::Row& r) {
    const double a = (double)r.m_eeg * r.m_eeg - (double)r.mee * r.mee;
    if (!(a > 0)) return -1.0;
    const double f = (scan::kMpi0 * scan::kMpi0 - (double)r.mee * r.mee) / a;
    const double m2 = (double)r.m4 * r.m4 + f * ((double)r.m5 * r.m5 - (double)r.m4 * r.m4);
    return m2 > 0 ? std::sqrt(m2) : -1.0;
}
inline bool Tagged(const scan::Row& r) { return r.g_idx > 0 && r.mee < kMeeMax && r.m_eeg >= kPi0Lo && r.m_eeg < kPi0Hi; }
inline scan::PeakCfg TagCfg() {
    scan::PeakCfg c;
    c.mean = kTagMean; c.sigma = kTagSigma; c.lo = kTagLo; c.hi = kTagHi; c.order = kTagOrder;
    c.winLo = c.mean - scan::kWinNsig * c.sigma; c.winHi = c.mean + scan::kWinNsig * c.sigma;
    return c;
}
inline scan::FillOpt TagOpt() {
    scan::FillOpt o;
    o.nbin = kTagNbin; o.lo = kTagHLo; o.hi = kTagHHi; o.reb = kTagReb;
    o.x = [](const scan::Row& r) { return M5c(r); };
    return o;
}

// ---- selections -------------------------------------------------------------------------------
inline bool Sel(int s, const scan::Row& r) { return PassSel(s, r, r.cor_ok == 1, r.dz_pi, r.mm4); }

// ---- F10 cell of a row (identical to wp0_fill.C: fine bin first, then its centre) ---------------
std::vector<double> gF10e;
inline int F10Cell(const scan::Row& r) {
    const int fi = FineBin(r.m4);
    if (fi < 0) return -1;
    const double c = kFineLo + (fi + 0.5) * kFineW;
    if (c < gF10e.front() || c >= gF10e.back()) return -1;
    const int x = (int)(std::upper_bound(gF10e.begin(), gF10e.end(), c) - gF10e.begin()) - 1;   // 0..69
    const int s = Slice(r.mpipi), e = Eclass(r.mee);
    if (s < 0 || e < 0) return -1;
    return (2 * s + e) * kNF10 + x;
}

// ---- spread-rule CB for a slice: CB(all slices) x f(slice | 50 MeV band) -----------------------
// bands of 50 MeV from 0.30 (WP0 grid); the engine's fine grid starts at 0.35 (aligned)
inline int BandOfX(double x) { return (int)std::floor((x - 0.30) / 0.05 + 1e-9); }
scan::Spectra SpreadCB(const scan::Spectra& all, const scan::Spectra& sl) {
    scan::Spectra s = sl;
    std::map<int, double> la, ls;
    for (int i = 0; i < sl.nbin; ++i) {
        const int b = BandOfX(sl.x(i));
        la[b] += all.n[1][i] + all.n[2][i];
        ls[b] += sl.n[1][i] + sl.n[2][i];
    }
    auto f = [&](double x) { const int b = BandOfX(x); return la[b] > 0 ? ls[b] / la[b] : 0.; };
    for (int i = 0; i < s.nbin; ++i) {
        const double fi = f(s.x(i));
        s.cb[i] = all.cb[i] * fi;
        s.cbVar[i] = all.cbVar[i] * fi * fi;
        s.sig[i] = s.n[0][i] - s.cb[i];
        s.sigVar[i] = s.n[0][i] + s.cbVar[i];
    }
    for (int j = 0; j < s.nw; ++j) {
        const double xc = s.lo + (j + 0.5) * s.reb * s.bw;
        const double fj = f(xc);
        s.cbW[j] = all.cbW[j] * fj;
        s.cbWVar[j] = all.cbWVar[j] * fj * fj;
    }
    return s;
}
inline double SumRange(const std::vector<double>& v, const scan::Spectra& sp, double a, double b) {
    double s = 0;
    for (int i = 0; i < sp.nbin; ++i) if (sp.x(i) >= a && sp.x(i) < b) s += v[i];
    return s;
}

}  // namespace wp3r

void wp3_rows_cnt() {
    using namespace alldefs;
    using namespace wp3r;
    gSystem->mkdir(kOutDir, kTRUE);
    gOut = fopen(Form("%s/wp3_rows.txt", kOutDir), "w");
    gSum = fopen(Form("%s/wp3_rows_summary.txt", kOutDir), "w");
    gF10e = F10Edges();
    P("wp3_rows.C  (WP3 stage A, row-level parts)  data: %s\n\n", scan::kCompact);

    scan::Data d;
    if (!scan::Load(d)) return;

    // ------------------------------------------------------------------------------------------
    // 1. event-group covariance of N over the F10 cells
    // ------------------------------------------------------------------------------------------
    const int nC = 10 * kNF10;
    const int selList[2] = {kSelS0, kSelPre};
    TFile* fo = TFile::Open(Form("%s/wp3_rows.root", kOutDir), "RECREATE");
    TFile* fh = TFile::Open(kHistFile);
    for (int q = 0; q < 2; ++q) {
        const int s = selList[q];
        TMatrixDSym S(nC);
        std::vector<double> cnt(nC, 0.);
        std::vector<int> cells;
        int cur = -1;
        long ngroups = 0, nmulti = 0;
        auto flush = [&]() {
            if (cells.empty()) return;
            ++ngroups;
            if (cells.size() > 1) ++nmulti;
            for (size_t a = 0; a < cells.size(); ++a)
                for (size_t b = 0; b < cells.size(); ++b) S(cells[a], cells[b]) += 1.0;
            cells.clear();
        };
        for (const auto& r : d.us()) {
            if (r.grp != cur) { flush(); cur = r.grp; }
            if (!Sel(s, r)) continue;
            const int c = F10Cell(r);
            if (c < 0) continue;
            cells.push_back(c);
            cnt[c] += 1;
        }
        flush();
        // check against hists_v1.root
        double maxd = 0, sumN = 0, sumVarDiag = 0, sumOff = 0;
        for (int h = 0; h < 10; ++h) {
            TH1* N = allio::Data(fh, SelTag(s), "all", "N", F10Name(h / 2, h % 2));
            for (int x = 0; x < kNF10; ++x) {
                const int c = h * kNF10 + x;
                maxd = std::max(maxd, std::fabs(N->GetBinContent(x + 1) - cnt[c]));
                maxd = std::max(maxd, std::fabs(N->GetBinError(x + 1) * N->GetBinError(x + 1) - S(c, c)));
                sumN += cnt[c];
            }
            delete N;
        }
        for (int a = 0; a < nC; ++a) { sumVarDiag += S(a, a); for (int b = 0; b < nC; ++b) if (a != b) sumOff += S(a, b); }
        P("1. event-group covariance %s: %.0f rows in F10, %ld groups (%ld with >1 row in F10); max |N - N_F10(hists_v1)| and |Sigma_aa - err^2| = %.3g\n",
          SelTag(s), sumN, ngroups, nmulti, maxd);
        P("   sum of the diagonal %.0f, sum of the off-diagonal %.0f: variance of the F10 total / N = %.4f\n", sumVarDiag, sumOff,
          (sumVarDiag + sumOff) / sumN);
        fo->cd();
        S.Write(Form("egcov_%s", SelTag(s)));
    }

    // ------------------------------------------------------------------------------------------
    // 2. A2 edge step in CR-K slices
    // ------------------------------------------------------------------------------------------
    const double eK = Edge(kEdgeLimit, 0.);
    P("\n2. A2 CR-K null test: pol2 + h x step at the FIXED edge E(%.4f, 0) = %.5f GeV/c2, smearing %.1f MeV, range edge -+ 45 MeV\n",
      kEdgeLimit, eK, 1e3 * scan::kEdgeSigma);
    P("   h = step height in rows per 2.5 MeV (positive = more rows below the edge). Earlier: none, 0.42-0.46, edge 0.5150: 36.5 +- 19.3\n");
    struct SlDef { double lo, hi; bool ctrl; };
    const std::vector<SlDef> sl = {{0.43, 0.46, false}, {0.46, 0.55, false}, {0.44, 0.46, false}, {0.43, 0.55, false}, {0.42, 0.46, true}};
    const int sels3[3] = {kSelS0, kSelPre, kSelNone};
    std::vector<scan::Sel> vs;
    for (int q = 0; q < 3; ++q) {
        const int s = sels3[q];
        vs.push_back([s](const scan::Row& r) { return Sel(s, r); });
        for (const auto& z : sl) {
            const float a = z.lo, b = z.hi;
            vs.push_back([s, a, b](const scan::Row& r) { return Sel(s, r) && r.mpipi >= a && r.mpipi < b; });
        }
    }
    const std::vector<scan::Spectra> sp = scan::FillMany(d, vs);
    const int nPer = 1 + (int)sl.size();
    struct EdgeRes { scan::StepFit own, spr; double nBelow, cbOwnBelow, cbSprBelow; };
    std::vector<EdgeRes> er(3 * sl.size());
    P("   %-5s %-9s %-6s | %-28s | %-28s | %s\n", "sel", "M(pipi)", "edge", "CB own: h +- err (h/err) chi2/ndf",
      "CB spread: h +- err (h/err) chi2/ndf", "N, CB own, CB spread for M4 < edge + 2 sigma");
    for (int q = 0; q < 3; ++q) {
        const scan::Spectra& all = sp[q * nPer];
        for (size_t k = 0; k < sl.size(); ++k) {
            const scan::Spectra& own = sp[q * nPer + 1 + k];
            const scan::Spectra spr = SpreadCB(all, own);
            scan::StepCfg c;
            c.mpLo = sl[k].lo; c.mpHi = sl[k].hi; c.useEdge = true; c.edge = sl[k].ctrl ? 0.5150 : eK;
            EdgeRes& R = er[q * sl.size() + k];
            R.own = scan::FitEdgeStep(own, c);
            R.spr = scan::FitEdgeStep(spr, c);
            const double lim = c.edge + 2 * scan::kEdgeSigma;
            R.nBelow = SumRange(own.n[0], own, 0.30, lim);
            R.cbOwnBelow = SumRange(own.cb, own, 0.30, lim);
            R.cbSprBelow = SumRange(spr.cb, spr, 0.30, lim);
            P("   %-5s %.2f-%.2f %.4f | %6.2f +- %5.2f (%5.2f) %5.1f/%-3d | %6.2f +- %5.2f (%5.2f) %5.1f/%-3d | %.0f, %.1f, %.1f%s\n",
              SelTag(sels3[q]), sl[k].lo, sl[k].hi, c.edge, R.own.h, R.own.eh, R.own.eh > 0 ? R.own.h / R.own.eh : 0., R.own.chi2,
              R.own.ndf, R.spr.h, R.spr.eh, R.spr.eh > 0 ? R.spr.h / R.spr.eh : 0., R.spr.chi2, R.spr.ndf, R.nBelow, R.cbOwnBelow,
              R.cbSprBelow, sl[k].ctrl ? "   (old control, edge by hand)" : "");
            if (gSum) fprintf(gSum, "A2row %s %.2f %.2f %.5f own %.4f %.4f %.3f %d spread %.4f %.4f %.3f %d nbelow %.0f\n", SelTag(sels3[q]), sl[k].lo, sl[k].hi,
                              c.edge, R.own.h, R.own.eh, R.own.chi2, R.own.ndf, R.spr.h, R.spr.eh, R.spr.chi2, R.spr.ndf, R.nBelow);
        }
        // check: spread-rule CB of S5 in CR-K = WP0 CB_F10_S5E1 + S5E2 in CR-K
        if (q < 2) {
            const scan::Spectra spr = SpreadCB(all, sp[q * nPer + 1 + 3]);   // 0.43-0.55
            const scan::Spectra& own = sp[q * nPer + 1 + 3];
            double wp0 = 0, n0 = 0;
            for (int e = 0; e < 2; ++e) {
                TH1* cb = allio::Data(fh, SelTag(sels3[q]), "all", "CB", F10Name(4, e));
                TH1* nn = allio::Data(fh, SelTag(sels3[q]), "all", "N", F10Name(4, e));
                for (int x = 1; x <= kNF10; ++x) if (cb->GetXaxis()->GetBinUpEdge(x) <= 0.5225 + 1e-9) { wp0 += cb->GetBinContent(x); n0 += nn->GetBinContent(x); }
                delete cb; delete nn;
            }
            P("   check %s: CB(spread) of 0.43-0.55 in M4 0.40-0.5225 = %.2f (N %.0f); WP0 CR-K CB of S5 (>= 0.43, all M(ee)) = %.2f (N %.0f)\n",
              SelTag(sels3[q]), SumRange(spr.cb, spr, 0.40, 0.5225), SumRange(own.n[0], own, 0.40, 0.5225), wp0, n0);
        }
    }

    // ------------------------------------------------------------------------------------------
    // 3. ECAL-tagged eta in CR-K
    // ------------------------------------------------------------------------------------------
    P("\n3. ECAL-tagged eta (M(e+e-gamma) %.3f-%.3f, M(ee) < %.2f): m5c peak content T (frozen %.4f / %.1f MeV, pol%d, %.2f-%.2f; P_main)\n",
      kPi0Lo, kPi0Hi, kMeeMax, kTagMean, 1e3 * kTagSigma, kTagOrder, kTagLo, kTagHi);
    std::vector<scan::Sel> vt;
    const char* tname[3] = {"CR-K (M(pipi) >= 0.43, M4 0.400-0.5225)", "bump region (M(pipi) < 0.43, M4 < 0.5225)", "all tagged rows"};
    for (int q = 0; q < 3; ++q) {
        const int s = sels3[q];
        vt.push_back([s](const scan::Row& r) { return Tagged(r) && Sel(s, r) && r.mpipi >= 0.43f && r.m4 >= 0.400f && r.m4 < 0.5225f; });
        vt.push_back([s](const scan::Row& r) { return Tagged(r) && Sel(s, r) && r.mpipi < 0.43f && r.m4 < 0.5225f; });
        vt.push_back([s](const scan::Row& r) { return Tagged(r) && Sel(s, r); });
    }
    const std::vector<scan::Spectra> tp = scan::FillMany(d, vt, TagOpt());
    std::vector<scan::PeakFit> tf(tp.size());
    for (int q = 0; q < 3; ++q)
        for (int k = 0; k < 3; ++k) {
            const scan::Spectra& s = tp[3 * q + k];
            tf[3 * q + k] = scan::FitPeakMain(s, TagCfg());
            const scan::PeakFit& f = tf[3 * q + k];
            P("   %-5s %-46s rows %6.0f (N++ %4.0f, N-- %4.0f): T = %8.1f +- %6.1f (T/err %5.2f), S/B %.3f, chi2/ndf %5.1f/%d\n", SelTag(sels3[q]),
              tname[k], s.nsel[0], s.nsel[1], s.nsel[2], f.S, f.eS, f.eS > 0 ? f.S / f.eS : 0., f.SoB, f.chi2, f.ndf);
            if (gSum) fprintf(gSum, "TAG %s %d %.0f %.3f %.3f %.4f %.3f %d\n", SelTag(sels3[q]), k, s.nsel[0], f.S, f.eS, f.SoB, f.chi2, f.ndf);
        }

    // ------------------------------------------------------------------------------------------
    // figures
    // ------------------------------------------------------------------------------------------
    bla::Style();
    gStyle->SetPadLeftMargin(0.16);
    {
        // edge-step fits in the CR-K slices, S0: 2 slices x 2 CB constructions
        TCanvas* c = new TCanvas("c_edge", "", 1400, 1000);
        c->Divide(2, 2, 0.002, 0.002);
        int pad = 0;
        for (int k = 0; k < 2; ++k)
            for (int v = 0; v < 2; ++v) {
                c->cd(++pad);
                const scan::Spectra& own = sp[0 * nPer + 1 + k];
                const scan::Spectra spr = SpreadCB(sp[0], own);
                const scan::Spectra& S = v == 0 ? spr : own;
                const scan::StepFit& f = v == 0 ? er[k].spr : er[k].own;
                TH1D* hn = new TH1D(Form("hn%d", pad), "", S.nbin, S.lo, S.hi);
                TH1D* hc = new TH1D(Form("hc%d", pad), "", S.nbin, S.lo, S.hi);
                TH1D* hs = new TH1D(Form("hs%d", pad), "", S.nbin, S.lo, S.hi);
                for (int i = 0; i < S.nbin; ++i) {
                    hn->SetBinContent(i + 1, S.n[0][i]); hn->SetBinError(i + 1, std::sqrt(S.n[0][i]));
                    hc->SetBinContent(i + 1, S.cb[i]); hc->SetBinError(i + 1, std::sqrt(S.cbVar[i]));
                    hs->SetBinContent(i + 1, S.sig[i]); hs->SetBinError(i + 1, std::sqrt(S.sigVar[i]));
                }
                hn->GetXaxis()->SetRangeUser(0.42, 0.56);
                hn->SetMinimum(-15); hn->SetMaximum(std::max(10., 1.35 * hn->GetMaximum()));
                hn->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
                hn->GetYaxis()->SetTitle("counts / 2.5 MeV/c^{2}");
                hn->SetMarkerStyle(20); hn->SetMarkerColor(kBlack); hn->SetLineColor(kBlack); hn->SetMarkerSize(0.9);
                hc->SetMarkerStyle(21); hc->SetMarkerColor(kRed + 1); hc->SetLineColor(kRed + 1); hc->SetMarkerSize(0.8);
                hs->SetMarkerStyle(20); hs->SetMarkerColor(kBlue + 1); hs->SetLineColor(kBlue + 1); hs->SetMarkerSize(0.9);
                hn->Draw("E1");
                hc->Draw("E1 same");
                hs->Draw("E1 same");
                if (f.ok) {
                    TGraph* g = new TGraph(), *g0 = new TGraph();
                    for (size_t i = 0; i < f.tx.size(); ++i) { g->SetPoint(i, f.tx[i], f.fy[i]); g0->SetPoint(i, f.tx[i], f.fy0[i]); }
                    g->SetLineColor(kGreen + 2); g->SetLineWidth(3); g->Draw("L");
                    g0->SetLineColor(kGreen + 2); g0->SetLineStyle(2); g0->SetLineWidth(2); g0->Draw("L");
                }
                TLine* l = new TLine(eK, -15, eK, hn->GetMaximum()); l->SetLineStyle(3); l->SetLineWidth(2); l->Draw();
                TLatex t; t.SetNDC(); t.SetTextSize(0.050);
                t.DrawLatex(0.175, 0.86, Form("S0, M(#pi^{+}#pi^{-}) %.2f-%.2f, CB %s", sl[k].lo, sl[k].hi, v == 0 ? "spread (nominal)" : "own 2#sqrt{N_{++}N_{--}}"));
                t.DrawLatex(0.175, 0.79, Form("h = %.2f #pm %.2f counts/2.5 MeV, #chi^{2}/ndf %.1f/%d", f.h, f.eh, f.chi2, f.ndf));
                if (pad == 1) {
                    TLegend* lg = new TLegend(0.19, 0.50, 0.55, 0.75);
                    lg->AddEntry(hn, "N", "pe"); lg->AddEntry(hc, "CB", "pe"); lg->AddEntry(hs, "N - CB", "pe");
                    lg->SetTextSize(0.045); lg->Draw();
                }
            }
        bla::SaveFig(c, Form("%s/wp3_crk_edge_cnt", kFigDir));
    }
    {
        // tagged eta: CR-K vs bump region, S0
        TCanvas* c = new TCanvas("c_tag", "", 1400, 560);
        c->Divide(2, 1, 0.002, 0.002);
        for (int k = 0; k < 2; ++k) {
            c->cd(k + 1);
            const scan::Spectra& S = tp[k];
            const scan::PeakFit& f = tf[k];
            TH1D* hs = new TH1D(Form("ts%d", k), "", S.nbin, S.lo, S.hi);
            for (int i = 0; i < S.nbin; ++i) { hs->SetBinContent(i + 1, S.sig[i]); hs->SetBinError(i + 1, std::sqrt(S.sigVar[i])); }
            hs->GetXaxis()->SetRangeUser(0.45, 0.66);
            hs->SetMinimum(k == 0 ? -4 : 0);
            hs->SetMaximum(k == 0 ? 12 : 1.3 * hs->GetMaximum());
            hs->GetXaxis()->SetTitle("m_{5c} = M(#pi^{+}#pi^{-}e^{+}e^{-}#gamma)_{#pi^{0}-constrained} [GeV/c^{2}]");
            hs->GetYaxis()->SetTitle("N - CB / 2.5 MeV/c^{2}");
            hs->SetMarkerStyle(20); hs->SetMarkerColor(kBlue + 1); hs->SetLineColor(kBlue + 1);
            hs->Draw("E1");
            if (f.ok) {
                TGraph* g = new TGraph(), *g0 = new TGraph();
                int n = 0;
                for (double m = kTagLo + 0.5 * S.bw; m < kTagHi; m += S.bw) { g->SetPoint(n, m, f.Total(m)); g0->SetPoint(n++, m, f.Bkg(m)); }
                g->SetLineColor(kGreen + 2); g->SetLineWidth(3); g->Draw("L");
                g0->SetLineColor(kGreen + 2); g0->SetLineStyle(2); g0->Draw("L");
            }
            TLatex t; t.SetNDC(); t.SetTextSize(0.048);
            t.DrawLatex(0.19, 0.86, k == 0 ? "S0, ECAL-tagged, CR-K (M(#pi^{+}#pi^{-}) #geq 0.43)" : "S0, ECAL-tagged, M(#pi^{+}#pi^{-}) < 0.43, M4 < 0.5225");
            t.DrawLatex(0.19, 0.79, Form("T = %.1f #pm %.1f", f.S, f.eS));
        }
        bla::SaveFig(c, Form("%s/wp3_crk_tag", kOutDir));   // _cnt copy: scratch only
    }
    fo->Close();
    fh->Close();
    fclose(gOut); gOut = nullptr;
    fclose(gSum); gSum = nullptr;
    printf("written %s/wp3_rows.{root,txt}\n", kOutDir);
}
