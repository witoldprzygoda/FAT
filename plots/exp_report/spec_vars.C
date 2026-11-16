// ============================================================================
// spec_vars.C  --  distributions of the kinematic variables of the current
//                  production: N, CB and N - CB for the selection pre, for all
//                  M(pi+pi-e+e-) and restricted to the peak and bump windows
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Production pippimepem_isBestAll_T4530 (FAT: event vertex z > -500 mm, PT3
// trigbit == 8192, all stored combinations), corrected momenta (cor_ok == 1).
//
// What it does
//   For ten variables x -- M(e+e-), M(pi+pi-), OA(e+e-), OA_lab((pi+pi-),(e+e-)),
//   the (pi+pi-)-(e+e-) angle in the rest frame forced to m_eta ("eta-rest
//   angle"), MM4 = missing mass of pi+pi-e+e-, pT and cms rapidity y of the
//   pi+pi-e+e- system, event vertex z, dz_pi = max(|z_pi+ - vz|, |z_pi- - vz|) --
//   it fills the unlike-sign rows N(x) and the CB(x), and N - CB, for
//     (a) pre, all M4;  (b) pre and M4 in the data peak window 0.53364-0.55416;
//     (c) pre and M4 in the data bump window 0.4550-0.5225 GeV/c2.
//   For the vertex z the base is pre WITHOUT vz > -160 mm (the cut is marked).
//   CB(x) (main, "M4-binned"): the CB is defined per M4 bin as everywhere in the
//   report, CB_j = 2 sqrt(N++_j N--_j) in 10 MeV M4 bins j (all x), interpolated
//   to 2.5 MeV and summed over the region with the coefficients c_j of ScanEngine
//   (CB_region = sum_j c_j CB_j). It is distributed in x by weighting every
//   like-sign row of M4 bin j with c_j sqrt(N--_j/N++_j) (e+e+ row) or
//   c_j sqrt(N++_j/N--_j) (e-e- row): the weights of bin j add up to c_j CB_j,
//   so sum_x CB(x) = CB_region exactly (the same number as in spec_counts.txt).
//   Check (table only): the bin-by-bin CB'(x) = 2 sqrt(N++(x) N--(x)) per x bin;
//   sum_x CB'(x) is smaller when the e+e+ and e-e- shapes in x differ.
//
// Selection pre = cor_ok && isBest >= 0 && vz > -160 mm (spec_common.h); the
// FAT production already has vertex z > -500 mm and PT3. Rows are particle
// combinations, not events; not normalised to luminosity.
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root
//          (scan::kCompact) through ../ScanEngine.h and spec_common.h
// Outputs: /home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/spec_var_<x>.pdf (+ .png), x = mee, mpipi,
//          oa, oalab, oaeta, mm4, pt4, y4, vz, dzpi
//          /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/spec_vars.txt
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/spec/spec_vars.root
//
// Run (from the repository root /home/przygoda/HADES/FAT/FAT; about 1 min, 4 GB RAM):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_spec", kTRUE)' 'plots/exp_report/spec_vars.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "spec_common.h"

namespace specvars {

using namespace spec;
using scan::Spectra;

struct VarDef {
    const char* tag;        // file tag
    const char* name;       // plain name (table)
    const char* axis;       // TLatex axis title
    const char* unit;       // bin unit for the y title
    scan::Var x;            // value
    double lo, hi; int nb;  // all M4
    double wlo, whi; int wnb;   // window panels
    bool logA, logW;        // log y: all-M4 panel, window panels
    bool legLeft;           // legend in the upper left of panel (a)
    bool noVz;              // base without vz > -160 (vertex z)
};

FILE* gOut = nullptr;
void P(const char* fmt, ...) {
    va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
    if (gOut) { va_start(a, fmt); vfprintf(gOut, fmt, a); va_end(a); }
}

// y range so that the data fill the lower 60 % of the frame
void Range(const TH1D* hn, const TH1D* hc, bool logy, double& y1, double& y2) {
    double mx = 0, mn = 1e30;
    for (int i = 1; i <= hn->GetNbinsX(); ++i) {
        mx = std::max(mx, hn->GetBinContent(i) + hn->GetBinError(i));
        if (hc->GetBinContent(i) > 0) mn = std::min(mn, hc->GetBinContent(i));
        if (hn->GetBinContent(i) > 0) mn = std::min(mn, hn->GetBinContent(i));
    }
    if (logy) {
        y1 = std::max(0.5, 0.5 * mn);
        const double l1 = std::log10(y1), l0 = std::log10(mx);
        y2 = std::pow(10., l1 + (l0 - l1) / 0.62);
    } else {
        y1 = 0;
        y2 = mx / 0.62;
    }
}

} // namespace specvars

void spec_vars() {
    using namespace specvars;
    TStopwatch sw;
    Style();
    gSystem->mkdir(kWork, kTRUE);
    scan::Data d;
    if (!scan::Load(d)) return;

    std::vector<VarDef> vars = {
        {"mee", "M(e+e-) [MeV/c2]", "M(e^{+}e^{-}) [MeV/c^{2}]", "5 MeV/c^{2}",
         [](const Row& r) { return 1000.0 * r.mee; }, 0, 600, 120, 0, 300, 60, true, true, false, false},
        {"mpipi", "M(pi+pi-) [GeV/c2]", "M(#pi^{+}#pi^{-}) [GeV/c^{2}]", "10 MeV/c^{2}",
         [](const Row& r) { return (double)r.mpipi; }, 0.25, 1.25, 100, 0.27, 0.55, 28, true, false, false, false},
        {"oa", "OA(e+e-) [deg]", "OA(e^{+}e^{-}) [deg]", "1 deg",
         [](const Row& r) { return (double)r.oa; }, 0, 90, 90, 0, 90, 45, true, true, false, false},
        {"oalab", "OA_lab((pi+pi-),(e+e-)) [deg]", "#alpha_{lab}(#pi^{+}#pi^{-}, e^{+}e^{-}) [deg]", "2 deg",
         [](const Row& r) { return (double)r.oa_lab; }, 0, 150, 75, 0, 150, 50, false, false, false, false},
        {"oaeta", "eta-rest angle ((pi+pi-),(e+e-)) [deg]", "#alpha_{#eta rest}(#pi^{+}#pi^{-}, e^{+}e^{-}) [deg]", "2 deg",
         [](const Row& r) { return (double)r.oa_eta; }, 0, 180, 90, 140, 180, 40, false, true, true, false},
        {"mm4", "MM(pi+pi-e+e-) [GeV/c2]", "MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]", "20 MeV/c^{2}",
         [](const Row& r) { return (double)r.mm4; }, 1.6, 3.2, 80, 2.0, 3.1, 44, false, false, true, false},
        {"pt4", "pT(pi+pi-e+e-) [GeV/c]", "p_{T}(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c]", "20 MeV/c",
         [](const Row& r) { return 1e-3 * r.pt4; }, 0, 1.2, 60, 0, 1.0, 40, false, false, false, false},
        {"y4", "y_cm(pi+pi-e+e-)", "y_{cm}(#pi^{+}#pi^{-}e^{+}e^{-})", "0.025",
         [](const Row& r) { return (double)r.y4; }, -1.0, 0.8, 72, -1.0, 0.8, 36, false, false, false, false},
        {"vz", "event vertex z [mm]", "event vertex v_{z} [mm]", "5 mm",
         [](const Row& r) { return (double)r.vz; }, -500, 100, 120, -300, 50, 70, true, true, false, true},
        {"dzpi", "dz_pi [mm]", "#Deltaz_{#pi} = max|z_{#pi} #minus v_{z}| [mm]", "2 mm",
         [](const Row& r) { return (double)r.dz_pi; }, 0, 200, 100, 0, 200, 50, true, true, false, false},
    };

    gOut = fopen(Form("%s/spec_vars.txt", kTabOut), "w");
    P("# spec_vars.txt -- kinematic distributions of the current production (written by plots/exp_report/spec_vars.C)\n");
    P("# data: %s ; FAT: vertex z > -500 mm, PT3 (trigbit 8192), all stored combinations\n", scan::kCompact);
    P("# selection pre = cor_ok && isBest >= 0 && vz > -160 mm (vertex z: base without vz > -160); regions in M4 [GeV/c2]:\n");
    P("#   all, peak window %.5f-%.5f, bump window %.4f-%.4f (data windows)\n", scan::kPeakLo, scan::kPeakHi, scan::kBumpLo, scan::kBumpHi);
    P("# CB(x): M4-binned CB (2 sqrt(N++ N--) per 10 MeV M4 bin, interpolated, k = 1) distributed in x with the like-sign rows\n");
    P("#   (weights c_j sqrt(N--_j/N++_j) and c_j sqrt(N++_j/N--_j)). Columns: N in the x range (rows), fraction of the region rows\n");
    P("#   outside the x range, sumCB = sum of CB(x) in the x range, sumCB_all = incl. under/overflow, CB_M4 = region CB of the\n");
    P("#   M4 method (spec_counts.txt), closure = sumCB_all / CB_M4 (must be 1), GM_x = sum over x bins of 2 sqrt(N++(x) N--(x))\n");
    P("#   of the region rows (bin-by-bin alternative), GM_x / 2 sqrt(sum N++ sum N--), CB/N in range, means of N, CB, N - CB.\n");
    P("# The means are taken INSIDE the x range of the histogram (rows outside it are ignored; <x>_N-CB from the bin centres).\n");
    P("#   Long tails outside the range (OA > 90 deg, dz_pi > 200 mm, eta-rest angle < 140 deg in the window panels) shift them.\n#\n");
    P("%-8s %-6s %9s %7s %10s %10s %10s %7s %10s %7s %6s %9s %9s %9s\n", "var", "region", "N", "out", "sumCB", "sumCB_all", "CB_M4",
      "closure", "GM_x", "GM/tot", "CB/N", "<x>_N", "<x>_CB", "<x>_N-CB");

    TFile* fout = TFile::Open(Form("%s/spec_vars.root", kWork), "RECREATE");
    const char* regTag[3] = {"all", "peak", "bump"};
    const double regLo[3] = {0.25, scan::kPeakLo, scan::kBumpLo}, regHi[3] = {3.25, scan::kPeakHi, scan::kBumpHi};

    // ---- M4 spectra of the two bases (2.5 MeV, CB in 10 MeV bins as in spec_m4.C) and the
    //      CB coefficients c_j of the three regions (same arithmetic as scan::FillMany)
    scan::Sel bases[2] = {[](const Row& r) { return Pre(r); }, [](const Row& r) { return PreNoVz(r); }};
    scan::FillOpt om; om.lo = 0.25; om.hi = 3.25; om.nbin = 1200; om.reb = 4;
    om.windows.push_back(std::make_pair(0.25, 3.25));
    const std::vector<Spectra> m4 = scan::FillMany(d, {bases[0], bases[1]}, om);
    std::vector<double> cj[2][3], wgt[2][3][2];     // [base][region] c_j ; like-sign row weight [q = pp, mm]
    for (int b = 0; b < 2; ++b)
        for (int g = 0; g < 3; ++g) {
            const Spectra& sp = m4[b];
            cj[b][g].assign(sp.nw, 0.);
            for (int i = 0; i < sp.nbin; ++i) {
                const double u = sp.Frac(i, regLo[g], regHi[g]);
                if (u <= 0) continue;
                int j0, j1; double f; sp.Interp(i, j0, j1, f);
                cj[b][g][j0] += u * (1 - f) / sp.reb; cj[b][g][j1] += u * f / sp.reb;
            }
            for (int q = 0; q < 2; ++q) {
                wgt[b][g][q].assign(sp.nw, 0.);
                for (int j = 0; j < sp.nw; ++j) {
                    const double own = sp.w[1 + q][j], oth = sp.w[2 - q][j];
                    wgt[b][g][q][j] = (own > 0 && oth > 0) ? cj[b][g][j] * std::sqrt(oth / own) : 0.;
                }
            }
        }
    const double cbReg[2][3] = {{m4[0].win[2].cb, m4[0].win[0].cb, m4[0].win[1].cb}, {m4[1].win[2].cb, m4[1].win[0].cb, m4[1].win[1].cb}};

    for (const VarDef& v : vars) {
        const int b = v.noVz ? 1 : 0;
        const scan::Sel& base = bases[b];
        TH1D *hN[3], *hC[3], *hS[3], *hP[3], *hM[3];
        for (int g = 0; g < 3; ++g) {
            const double lo = g == 0 ? v.lo : v.wlo, hi = g == 0 ? v.hi : v.whi; const int nb = g == 0 ? v.nb : v.wnb;
            hN[g] = new TH1D(Form("%s_%s_us", v.tag, regTag[g]), "", nb, lo, hi);
            hC[g] = new TH1D(Form("%s_%s_cb", v.tag, regTag[g]), "", nb, lo, hi);
            hP[g] = new TH1D(Form("%s_%s_pp", v.tag, regTag[g]), "", nb, lo, hi);
            hM[g] = new TH1D(Form("%s_%s_mm", v.tag, regTag[g]), "", nb, lo, hi);
            for (TH1D* h : {hN[g], hC[g], hP[g], hM[g]}) { h->SetDirectory(nullptr); h->Sumw2(); }
        }
        double nSel[3] = {0, 0, 0};
        for (int k = 0; k < 3; ++k)
            for (const Row& r : d.s[k]) {
                if (!r.cor_ok || !base(r)) continue;
                const double m = r.m4;
                const double x = v.x(r);
                const int j = (int)((m - 0.25) / 0.01);
                for (int g = 0; g < 3; ++g) {
                    const bool inReg = (m >= regLo[g] && m < regHi[g]);
                    if (k == 0) {
                        if (!inReg) continue;
                        nSel[g] += 1;
                        hN[g]->Fill(x);
                    } else {
                        if (inReg) (k == 1 ? hP[g] : hM[g])->Fill(x);
                        if (j >= 0 && j < (int)wgt[b][g][k - 1].size() && wgt[b][g][k - 1][j] > 0) hC[g]->Fill(x, wgt[b][g][k - 1][j]);
                    }
                }
            }
        TCanvas* c = new TCanvas(Form("c_%s", v.tag), "", 2100, 700);
        c->Divide(3, 1, 0.001, 0.001);
        for (int g = 0; g < 3; ++g) {
            hS[g] = (TH1D*)hN[g]->Clone(Form("%s_%s_sig", v.tag, regTag[g]));
            hS[g]->SetDirectory(nullptr);
            hS[g]->Add(hC[g], -1.0);
            TH1D* hn = hN[g]; TH1D* hc = hC[g]; TH1D* hs = hS[g];
            Deco(hn, kBlack, 20, 0.7); Deco(hc, kRed + 1, 20, 0.7); Deco(hs, kBlue + 1, 20, 0.7);
            fout->cd(); hn->Write(); hc->Write(); hs->Write(); hP[g]->Write(); hM[g]->Write();
            // table numbers (in the histogram range)
            const int nb = hn->GetNbinsX();
            const double sN = hn->Integral(1, nb), sC = hc->Integral(1, nb), sCall = hc->Integral(0, nb + 1);
            const double sP = hP[g]->Integral(1, nb), sM = hM[g]->Integral(1, nb);
            double sGM = 0;
            for (int i = 1; i <= nb; ++i) sGM += 2 * std::sqrt(hP[g]->GetBinContent(i) * hM[g]->GetBinContent(i));
            const double cbt = 2 * std::sqrt(sP * sM);
            P("%-8s %-6s %9.0f %7.4f %10.1f %10.1f %10.1f %7.4f %10.1f %7.4f %6.3f %9.4g %9.4g %9.4g\n", v.tag, regTag[g], sN,
              nSel[g] > 0 ? 1 - sN / nSel[g] : 0., sC, sCall, cbReg[b][g], sCall / cbReg[b][g], sGM, cbt > 0 ? sGM / cbt : 0.,
              sN > 0 ? sC / sN : 0., hn->GetMean(), hc->GetMean(), hs->GetMean());
            // figure
            c->cd(g + 1);
            gPad->SetRightMargin(0.05);   // verifier: room for the last x tick label (was clipped at 0.03)
            const bool logy = (g == 0) ? v.logA : v.logW;
            if (logy) gPad->SetLogy();
            // linear panels with more than 2000 rows per bin: axis in units of 10^3 rows
            double y1, y2; Range(hn, hc, logy, y1, y2);
            const bool kilo = !logy && y2 > 2000.;
            if (kilo) {
                for (TH1D* h : {hn, hc, hs})
                    for (int i = 0; i <= h->GetNbinsX() + 1; ++i) { h->SetBinContent(i, 1e-3 * h->GetBinContent(i)); h->SetBinError(i, 1e-3 * h->GetBinError(i)); }
                y2 *= 1e-3;
            }
            TH1D* fr = new TH1D(Form("fr_%s_%d", v.tag, g), "", 100, hn->GetXaxis()->GetXmin(), hn->GetXaxis()->GetXmax());
            fr->SetDirectory(nullptr); fr->SetMinimum(y1); fr->SetMaximum(y2);
            fr->GetXaxis()->SetTitle(Fx(v.axis).c_str());
            fr->GetYaxis()->SetTitle(Form("%s / %s", kilo ? "10^{3} rows" : "rows", v.unit));
            fr->GetXaxis()->SetNdivisions(506);
            fr->Draw("axis");
            if (v.noVz) VLine(-160, y1, logy ? y1 * std::pow(y2 / y1, 0.6) : 0.6 * y2, kGray + 2, 2);
            hc->Draw("e same"); hn->Draw("e same"); hs->Draw("e same");
            gPad->RedrawAxis();
            const char* bt = v.noVz ? "no v_{z} cut" : "pre";
            if (g == 0) Txt(0.16, 0.935, Form("(a) %s, all M_{4}", bt), 0.050);
            if (g == 1) Txt(0.16, 0.935, Form("(b) %s, peak window", bt), 0.050);
            if (g == 2) Txt(0.16, 0.935, Form("(c) %s, bump window", bt), 0.050);
            if (g == 0) {
                TLegend* lg = v.legLeft ? Leg(0.20, 0.66, 0.60, 0.89, 0.046) : Leg(0.54, 0.66, 0.94, 0.89, 0.046);
                LE(lg, hn, "N (e^{+}e^{-} rows)", "pe");
                LE(lg, hc, "CB (M_{4}-binned)", "pe");
                LE(lg, hs, "N #minus CB", "pe");
                lg->Draw();
            } else {
                const double a = regLo[g], bb = regHi[g];
                const double xx = v.legLeft ? 0.20 : 0.93; const int al = v.legLeft ? 11 : 31;
                Txt(xx, 0.84, Form("%.4f #leq M_{4} < %.4f GeV/c^{2}", a, bb), 0.044, al);
                Txt(xx, 0.77, Form("CB/N = %.3f", sN > 0 ? sC / sN : 0.), 0.044, al);
            }
        }
        Save(c, Form("spec_var_%s", v.tag));
    }
    fout->Close();
    fclose(gOut); gOut = nullptr;
    printf("spec_vars done in %.1f s\n", sw.RealTime());
}
