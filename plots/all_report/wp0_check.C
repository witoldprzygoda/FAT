// ============================================================================
// wp0_check.C  --  re-reads hists_v1.root (WP0): acceptance values, closure,
// replica validation, definition figure, inventory tables (prefix wp0)
// ============================================================================
// Outputs
//   docs/latex/figs_all/wp0_definitions.{pdf,png}   layout F10, windows, CRs, edge E(m,q) on the simulated bump
//   docs/latex/figs_all/wp0_cb_errors.{pdf,png}     CB error per F10 cell: replica RMS vs linear propagation (data S0)
//   docs/latex/figs_all/wp0_sim_mcerr.{pdf,png}     sim MC error per F10 cell: replica RMS vs event-group error (S0__nom)
//   docs/latex/tables_all/wp0_check.txt             numbers read back from the file
//   docs/latex/tables_all/wp0_sim_totals.{txt,tex}  sim cross sections per directory (tot, B, P)
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp0", kTRUE)' 'plots/all_report/wp0_check.C+O'
// The sealed file is never opened here.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TBox.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"

#include "defs.h"
#include "hists_io.h"
#include "s5_blind_ana_fit.h"

namespace wp0c {
using namespace alldefs;
const char* const kFig = "/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all";
const char* const kTab = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all";

double RepRms(const TH1* rep, int bin, int nrep, double* mean = nullptr) {
    double m = 0, m2 = 0;
    for (int r = 1; r <= nrep; ++r) { const double x = rep->GetBinContent(bin, r); m += x; m2 += x * x; }
    m /= nrep;
    if (mean) *mean = m;
    return std::sqrt(std::max(0., m2 / nrep - m * m));
}
}  // namespace wp0c

void wp0_check() {
    using namespace wp0c;
    bla::Style();
    gStyle->SetPalette(kGreyScale);
    TColor::InvertPalette();
    TFile* f = TFile::Open(kHistFile);
    if (!f || f->IsZombie()) { printf("wp0_check: cannot open %s\n", kHistFile); return; }
    FILE* o = fopen(Form("%s/wp0_check.txt", kTab), "w");
    fprintf(o, "# wp0_check.txt -- plots/all_report/wp0_check.C, numbers read back from %s\n", kHistFile);

    // ---------------------------------------------------------------- acceptance values from the file
    {
        TH1* t = allio::Sim(f, "pre", "nom", kT12, "tot");
        fprintf(o, "A1  sim pre__nom T12 all M4 = %.4f +- %.4f nb (ref 701.1)\n", t->GetBinContent(1), t->GetBinError(1));
        TH1* w = allio::Data(f, "S0", "all", "N", "win");
        TH1* c = allio::Data(f, "S0", "all", "CB", "win");
        TH2* cr = (TH2*)allio::Data(f, "S0", "all", "CBrep", "win");
        const double rms = RepRms(cr, kWP + 1, kNRepCB);
        fprintf(o, "A2  data S0 peak window N = %.0f (ref 55947), CB = %.1f +- %.1f (replicas) / %.1f (linear), CB/N = %.4f (ref 0.059)\n",
                allio::Win(w, kWP), allio::Win(c, kWP), rms, allio::WinErr(c, kWP), allio::Win(c, kWP) / allio::Win(w, kWP));
        fprintf(o, "    data S0 X17 counting window X17c (M4 0.5312-0.5566, M(ee) 15.5-24.5 MeV) N = %.0f, CB = %.1f +- %.1f (replicas)\n",
                allio::Win(w, kWX17c), allio::Win(c, kWX17c), RepRms(cr, kWX17c + 1, kNRepCB));
        fprintf(o, "    data S0 window X17P (P and M(ee) 15.5-24.5 MeV) N = %.0f, CB = %.1f +- %.1f (replicas) (D19 quotes 5258 and CB 112; "
                   "res/smp notes 101.7)\n", allio::Win(w, kWX17P), allio::Win(c, kWX17P), RepRms(cr, kWX17P + 1, kNRepCB));
        TH1* wn = allio::Data(f, "none", "all", "N", "win");
        TH1* cn = allio::Data(f, "none", "all", "CB", "win");
        TH2* crn = (TH2*)allio::Data(f, "none", "all", "CBrep", "win");
        const double e = std::sqrt(std::pow(allio::WinErr(wn, kWB), 2) + std::pow(RepRms(crn, kWB + 1, kNRepCB), 2));
        fprintf(o, "A3  data none bump window W = N - CB = %.1f +- %.1f (ref 148091 +- 469); N = %.0f, CB = %.1f, CB/N = %.4f\n",
                allio::Win(wn, kWB) - allio::Win(cn, kWB), e, allio::Win(wn, kWB), allio::Win(cn, kWB), allio::Win(cn, kWB) / allio::Win(wn, kWB));
        delete t; delete w; delete c; delete cr; delete wn; delete cn; delete crn;
    }
    // ---------------------------------------------------------------- CB spread rule vs per-cell geometric mean (S1 information)
    {
        fprintf(o, "\n# CB of the spread rule vs 2 sqrt(N++ N--) of the same cell (summed over the cell), data S0 all\n");
        fprintf(o, "# %-10s %8s %8s %8s %12s %12s %8s\n", "cell", "N++", "N--", "N--/N++", "CB(spread)", "2sqrt(++--)", "ratio");
        for (int e = 0; e < 2; ++e)
            for (int s = 0; s < kNSlice; ++s) {
                TH1* pp = allio::Data(f, "S0", "all", "Npp", allio::F10(s, e));
                TH1* mm = allio::Data(f, "S0", "all", "Nmm", allio::F10(s, e));
                TH1* cb = allio::Data(f, "S0", "all", "CB", allio::F10(s, e));
                const double a1 = pp->Integral(1, kNF10), a2 = mm->Integral(1, kNF10), c = cb->Integral(1, kNF10), g = 2 * std::sqrt(a1 * a2);
                fprintf(o, "  F10_S%dE%d  %8.0f %8.0f %8.2f %12.1f %12.1f %8.3f\n", s + 1, e + 1, a1, a2, a2 / std::max(1., a1), c, g, g / c);
                delete pp; delete mm; delete cb;
            }
        TH1* pp = allio::Data(f, "S0", "all", "Npp", "win");
        TH1* mm = allio::Data(f, "S0", "all", "Nmm", "win");
        TH1* cb = allio::Data(f, "S0", "all", "CB", "win");
        for (int w : {kWB, kWP, kWX17c, kWX17P, kWCRU, kWCRK}) {
            const double a1 = allio::Win(pp, w), a2 = allio::Win(mm, w), c = allio::Win(cb, w), g = 2 * std::sqrt(a1 * a2);
            fprintf(o, "  win %-6s %8.0f %8.0f %8.2f %12.1f %12.1f %8.3f\n", kWin[w].tag, a1, a2, a2 / std::max(1., a1), c, g, g / std::max(1e-9, c));
        }
        delete pp; delete mm; delete cb;
    }
    // ---------------------------------------------------------------- closure from the file (F10 vs fine, all directories)
    {
        double worst = 0;
        std::string wname;
        int ntest = 0;
        const char* dsel[11] = {"S0", "pre", "dz40", "PEAK", "vz150", "vz170", "none", "S0_best", "pre_best", "S0_corfb", "pre_corfb"};
        for (const char* ds : dsel)
            for (int sb = 0; sb < kNSub; ++sb)
                for (const char* wh : {"N", "Npp", "Nmm", "CB"}) {
                    double s10 = 0;
                    for (int s = 0; s < kNSlice; ++s) for (int e = 0; e < 2; ++e) {
                        TH1* h = allio::Data(f, ds, SubTag(sb), wh, allio::F10(s, e));
                        s10 += h->Integral(1, kNF10); delete h;
                    }
                    TH1* fi = allio::Data(f, ds, SubTag(sb), wh, "fine_M4");
                    const double ref = fi->Integral(41, 160);
                    delete fi;
                    const double d = std::fabs(s10 - ref) / std::max(1., std::fabs(ref));
                    ++ntest;
                    if (d > worst) { worst = d; wname = Form("data %s %s %s", ds, SubTag(sb), wh); }
                }
        TDirectory* sd = f->GetDirectory("sim");
        for (TObject* k : *sd->GetListOfKeys()) {
            const std::string dn = k->GetName();
            const size_t p = dn.find("__");
            const std::string sel = dn.substr(0, p), var = dn.substr(p + 2);
            for (int comp = 0; comp < kNComp; ++comp) {
                double s10 = 0;
                for (int s = 0; s < kNSlice; ++s) for (int e = 0; e < 2; ++e) {
                    TH1* h = allio::Sim(f, sel.c_str(), var.c_str(), comp, allio::F10(s, e));
                    s10 += h->Integral(1, kNF10); delete h;
                }
                TH1* fi = allio::Sim(f, sel.c_str(), var.c_str(), comp, "fine_M4");
                const double ref = fi->Integral(41, 160);
                delete fi;
                const double d = std::fabs(s10 - ref) / std::max(1e-12, std::fabs(ref));
                ++ntest;
                if (ref != 0 && d > worst) { worst = d; wname = "sim " + dn + " " + CompTag(comp); }
            }
        }
        fprintf(o, "A4  closure read back from the file: sum of the 10 F10 histograms vs fine_M4 0.40-0.70, %d (directory, quantity) pairs: "
                   "max relative deviation %.2e (%s)\n", ntest, worst, wname.c_str());
    }
    // ---------------------------------------------------------------- sim totals per directory (table)
    {
        FILE* t = fopen(Form("%s/wp0_sim_totals.txt", kTab), "w");
        FILE* x = fopen(Form("%s/wp0_sim_totals.tex", kTab), "w");
        fprintf(t, "# wp0_sim_totals.txt -- plots/all_report/wp0_check.C; simulation, sigma [nb] (sum w / L_sim, L_sim = 221.05 nb^-1), "
                   "event-group errors; tot = all M4; B, P = data-scale windows on M4' (mass transfer of the variant)\n");
        fprintf(t, "%-16s %22s %22s %22s %22s %22s %22s\n", "directory", "T12 tot", "C2+C3 in B", "C4 in P", "C5 in P", "G_eta in P", "C6 (per BR) in X17c");
        fprintf(x, "%% wp0_sim_totals.tex -- plots/all_report/wp0_check.C (WP0); sigma [nb], event-group errors\n");
        fprintf(x, "\\begin{tabular}{lrrrrr}\n\\toprule\ndirectory & $T_{12}$ all $M_4$ & $C_2{+}C_3$ in B & $C_4$ in P & $C_5$ in P & $G_\\eta$ in P \\\\\n\\midrule\n");
        TDirectory* sd = f->GetDirectory("sim");
        for (TObject* k : *sd->GetListOfKeys()) {
            const std::string dn = k->GetName();
            const size_t p = dn.find("__");
            const std::string sel = dn.substr(0, p), var = dn.substr(p + 2);
            TH1* tt = allio::Sim(f, sel.c_str(), var.c_str(), kT12, "tot");
            TH1* w2 = allio::Sim(f, sel.c_str(), var.c_str(), kC2, "win");
            TH1* w3 = allio::Sim(f, sel.c_str(), var.c_str(), kC3, "win");
            TH1* w4 = allio::Sim(f, sel.c_str(), var.c_str(), kC4, "win");
            TH1* w5 = allio::Sim(f, sel.c_str(), var.c_str(), kC5, "win");
            TH1* we = allio::Sim(f, sel.c_str(), var.c_str(), kGeta, "win");
            TH1* w6 = allio::Sim(f, sel.c_str(), var.c_str(), kC6, "win");
            const double b23 = allio::Win(w2, kWB) + allio::Win(w3, kWB);
            const double e23 = std::hypot(allio::WinErr(w2, kWB), allio::WinErr(w3, kWB));
            fprintf(t, "%-16s %10.3f +- %8.3f %10.4f +- %8.4f %10.4f +- %8.4f %10.4f +- %8.4f %10.4f +- %8.4f %10.2f +- %8.2f\n", dn.c_str(),
                    tt->GetBinContent(1), tt->GetBinError(1), b23, e23, allio::Win(w4, kWP), allio::WinErr(w4, kWP), allio::Win(w5, kWP),
                    allio::WinErr(w5, kWP), allio::Win(we, kWP), allio::WinErr(we, kWP), allio::Win(w6, kWX17c), allio::WinErr(w6, kWX17c));
            std::string tag = dn;
            for (size_t q = 0; (q = tag.find('_', q)) != std::string::npos; q += 2) tag.replace(q, 1, "\\_");
            fprintf(x, "\\texttt{%s} & %.1f(%.0f) & %.3f(%.0f) & %.3f(%.0f) & %.4f(%.0f) & %.3f(%.0f) \\\\\n", tag.c_str(), tt->GetBinContent(1),
                    tt->GetBinError(1) * 10, b23, e23 * 1000, allio::Win(w4, kWP), allio::WinErr(w4, kWP) * 1000, allio::Win(w5, kWP),
                    allio::WinErr(w5, kWP) * 1e4, allio::Win(we, kWP), allio::WinErr(we, kWP) * 1000);
            delete tt; delete w2; delete w3; delete w4; delete w5; delete we; delete w6;
        }
        fprintf(x, "\\bottomrule\n\\end{tabular}\n");
        fclose(t); fclose(x);
    }
    // ---------------------------------------------------------------- figure 1: definitions on the simulated bump
    {
        TCanvas* c = new TCanvas("c1", "", 1000, 800);
        c->SetRightMargin(0.15);
        c->SetLeftMargin(0.13);
        TH2* m2 = (TH2*)allio::Sim(f, "S0", "nom", kC2, "B_map");
        TH2* m3 = (TH2*)allio::Sim(f, "S0", "nom", kC3, "B_map");
        m2->Add(m3);
        m2->Scale(1. / 0.01 / 0.01 * 1e-3);   // nb per bin -> ub / (GeV/c2)^2
        m2->GetXaxis()->SetRangeUser(0.40, 0.70);
        m2->GetYaxis()->SetRangeUser(0.27, 0.55);
        m2->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]");
        m2->GetYaxis()->SetTitle("M(#pi^{+}#pi^{-}) [GeV/c^{2}]");
        m2->GetZaxis()->SetTitle("d^{2}#sigma/dM_{4}dM_{#pi#pi} [#mub/(GeV/c^{2})^{2}]");
        m2->GetZaxis()->SetTitleOffset(1.3);
        m2->Draw("COLZ");
        // windows B and P (vertical bands), CR-U, CR-K, CR-L boxes
        auto box = [&](double x1, double y1, double x2, double y2, int col, int style) {
            TBox* b = new TBox(x1, y1, x2, y2);
            b->SetFillStyle(0); b->SetLineColor(col); b->SetLineWidth(3); b->SetLineStyle(style); b->Draw();
        };
        box(kWin[kWB].lo, 0.27, kWin[kWB].hi, 0.55, kBlue + 1, 1);
        box(kWin[kWP].lo, 0.27, kWin[kWP].hi, 0.55, kRed + 1, 1);
        box(kWin[kWCRU].lo, 0.27, kWin[kWCRU].hi, 0.55, kGreen + 2, 1);
        box(kWin[kWCRK].lo, 0.43, kWin[kWCRK].hi, 0.55, kMagenta + 1, 1);
        box(kWin[kWCRL].lo, 0.27, kWin[kWCRL].hi, 0.55, kOrange + 7, 2);
        for (int s = 0; s < kNSlice - 1; ++s) {
            TLine* l = new TLine(0.40, kSliceEdge[s], 0.70, kSliceEdge[s]);
            l->SetLineStyle(3); l->SetLineWidth(2); l->Draw();
        }
        // edge E(m, q) for q = 0 and 60 MeV
        for (int k = 0; k < 2; ++k) {
            TGraph* g = new TGraph();
            for (int i = 0; i <= 200; ++i) {
                const double m = 0.27 + (kEdgeLimit - 0.27) * i / 200.;
                const double e = Edge(m, k == 0 ? 0. : 0.060);
                if (e > 0) g->SetPoint(g->GetN(), e, m);
            }
            g->SetLineColor(k == 0 ? kBlack : kGray + 2); g->SetLineWidth(3); g->SetLineStyle(k == 0 ? 1 : 7);
            g->Draw("L");
        }
        TLatex tx; tx.SetTextSize(0.033);
        tx.DrawLatex(0.480, 0.535, "#color[602]{B}");
        tx.DrawLatex(0.537, 0.535, "#color[634]{P}");
        tx.DrawLatex(0.60, 0.535, "#color[418]{CR-U}");
        tx.DrawLatex(0.44, 0.45, "#color[617]{CR-K (S5)}");
        for (int s = 0; s < kNSlice; ++s) {
            const double y = s == 0 ? 0.30 : s == kNSlice - 1 ? 0.505 : 0.5 * (kSliceEdge[s - 1] + kSliceEdge[s]);
            tx.DrawLatex(0.680, y - 0.004, Form("S%d", s + 1));
        }
        tx.SetNDC(); tx.SetTextSize(0.030);
        tx.DrawLatex(0.13, 0.95, "simulation: bump C_{2}+C_{3}, S0, w#timesV3; M_{4}' = 0.99615 M_{4}, M_{#pi#pi}' = 0.9952 M_{#pi#pi}");
        TLegend* lg = new TLegend(0.54, 0.66, 0.80, 0.78);
        lg->SetFillStyle(1001);
        TGraph* g0 = new TGraph(); g0->SetLineWidth(3);
        TGraph* g1 = new TGraph(); g1->SetLineWidth(3); g1->SetLineStyle(7); g1->SetLineColor(kGray + 2);
        lg->AddEntry(g0, "edge E(m, 0)", "l");
        lg->AddEntry(g1, "edge E(m, 60 MeV/c^{2})", "l");
        lg->SetTextSize(0.032);
        lg->Draw();
        bla::SaveFig(c, Form("%s/wp0_definitions", kFig));
        delete c;
    }
    // ---------------------------------------------------------------- figure 2: CB error per F10 cell
    {
        TCanvas* c = new TCanvas("c2", "", 1600, 700);
        c->Divide(5, 2, 0.002, 0.002);
        fprintf(o, "\n# CB error per F10 histogram (data S0, all): sum over the 70 bins of CB; median over bins of (replica RMS / linear error)\n");
        for (int e = 0; e < 2; ++e)
            for (int s = 0; s < kNSlice; ++s) {
                c->cd(1 + s + 5 * e);
                gPad->SetLeftMargin(0.21); gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.04); gPad->SetRightMargin(0.03);
                TH1* cb = allio::Data(f, "S0", "all", "CB", allio::F10(s, e));
                TH2* rp = (TH2*)allio::Data(f, "S0", "all", "CBrep", allio::F10(s, e));
                TH1D* hl = new TH1D(Form("hl%d%d", s, e), "", kNF10, cb->GetXaxis()->GetXbins()->GetArray());
                TH1D* hr = (TH1D*)hl->Clone(Form("hr%d%d", s, e));
                std::vector<double> ratio;
                for (int b = 1; b <= kNF10; ++b) {
                    const double v = cb->GetBinContent(b);
                    if (v <= 0) continue;
                    const double r = RepRms(rp, b, kNRepCB);
                    hl->SetBinContent(b, cb->GetBinError(b) / v);
                    hr->SetBinContent(b, r / v);
                    ratio.push_back(r / std::max(1e-12, cb->GetBinError(b)));
                }
                std::sort(ratio.begin(), ratio.end());
                fprintf(o, "S%dE%d  CB %10.1f  median(rep/lin) %.3f\n", s + 1, e + 1, cb->Integral(1, kNF10), ratio.empty() ? 0. : ratio[ratio.size() / 2]);
                hr->SetLineColor(kRed + 1); hl->SetLineColor(kBlack); hl->SetLineStyle(2);
                hr->SetMaximum(0.55); hr->SetMinimum(0.);
                hr->GetXaxis()->SetTitle("M_{4} [GeV/c^{2}]");
                hr->GetYaxis()->SetTitle("#delta CB / CB");
                hr->GetXaxis()->SetLabelSize(0.06); hr->GetYaxis()->SetLabelSize(0.06);
                hr->GetXaxis()->SetTitleSize(0.065); hr->GetYaxis()->SetTitleSize(0.065);
                hr->GetXaxis()->SetNdivisions(505);
                hr->GetYaxis()->SetTitleOffset(1.45);
                hr->Draw("HIST");
                hl->Draw("HIST SAME");
                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.075);
                tx.DrawLatex(0.26, 0.87, Form("S%d E%d", s + 1, e + 1));
                if (s == 0 && e == 0) {
                    TLegend* lg = new TLegend(0.30, 0.58, 0.99, 0.80);
                    lg->AddEntry(hr, "200 replicas (RMS)", "l");
                    lg->AddEntry(hl, "linear, 10 MeV cells", "l");
                    lg->SetTextSize(0.06);
                    lg->Draw();
                }
                delete cb; delete rp;
            }
        bla::SaveFig(c, Form("%s/wp0_cb_errors", kFig));
        delete c;
    }
    // ---------------------------------------------------------------- figure 3: simulation MC error per F10 cell
    {
        TCanvas* c = new TCanvas("c3", "", 1600, 700);
        c->Divide(5, 2, 0.002, 0.002);
        fprintf(o, "\n# sim S0__nom: T12 and C5 per F10 histogram: sum [nb], N_eff, median over bins of (replica RMS / event-group error), "
                   "half A / half B\n");
        for (int e = 0; e < 2; ++e)
            for (int s = 0; s < kNSlice; ++s) {
                c->cd(1 + s + 5 * e);
                gPad->SetLeftMargin(0.21); gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.04); gPad->SetRightMargin(0.03);
                TH1* h = allio::Sim(f, "S0", "nom", kT12, allio::F10(s, e));
                TH2* rp = (TH2*)allio::Sim(f, "S0", "nom", kT12, allio::F10(s, e), "_rep");
                TH1* hA = allio::Sim(f, "S0", "nom", kT12, allio::F10(s, e), "_hA");
                TH1* hB = allio::Sim(f, "S0", "nom", kT12, allio::F10(s, e), "_hB");
                TH1* h5 = allio::Sim(f, "S0", "nom", kC5, allio::F10(s, e));
                TH1D* he = new TH1D(Form("he%d%d", s, e), "", kNF10, h->GetXaxis()->GetXbins()->GetArray());
                TH1D* hr = (TH1D*)he->Clone(Form("hq%d%d", s, e));
                std::vector<double> ratio;
                for (int b = 1; b <= kNF10; ++b) {
                    const double v = h->GetBinContent(b);
                    if (v <= 0) continue;
                    const double r = RepRms(rp, b, kNRepSim);
                    he->SetBinContent(b, h->GetBinError(b) / v);
                    hr->SetBinContent(b, r / v);
                    ratio.push_back(r / std::max(1e-30, h->GetBinError(b)));
                }
                std::sort(ratio.begin(), ratio.end());
                double eT = 0, e5 = 0;
                const double sT = h->IntegralAndError(1, kNF10, eT), s5 = h5->IntegralAndError(1, kNF10, e5);
                fprintf(o, "S%dE%d  T12 %9.4f nb (N_eff %8.0f)  C5 %8.5f nb (N_eff %8.0f)  median(rep/err) %.3f  hA/hB %.4f\n", s + 1, e + 1, sT,
                        eT > 0 ? sT * sT / eT / eT : 0., s5, e5 > 0 ? s5 * s5 / e5 / e5 : 0., ratio.empty() ? 0. : ratio[ratio.size() / 2],
                        hA->Integral(1, kNF10) / hB->Integral(1, kNF10));
                hr->SetLineColor(kRed + 1); he->SetLineColor(kBlack); he->SetLineStyle(2);
                hr->SetMaximum(0.55); hr->SetMinimum(0.);
                hr->GetXaxis()->SetTitle("M_{4}' [GeV/c^{2}]");
                hr->GetYaxis()->SetTitle("#delta T_{12} / T_{12} (MC)");
                hr->GetXaxis()->SetLabelSize(0.06); hr->GetYaxis()->SetLabelSize(0.06);
                hr->GetXaxis()->SetTitleSize(0.065); hr->GetYaxis()->SetTitleSize(0.065);
                hr->GetXaxis()->SetNdivisions(505);
                hr->GetYaxis()->SetTitleOffset(1.45);
                hr->Draw("HIST");
                he->Draw("HIST SAME");
                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.075);
                tx.DrawLatex(0.26, 0.87, Form("S%d E%d", s + 1, e + 1));
                if (s == 0 && e == 0) {
                    TLegend* lg = new TLegend(0.30, 0.58, 0.99, 0.80);
                    lg->AddEntry(hr, "100 replicas (RMS)", "l");
                    lg->AddEntry(he, "event-group error", "l");
                    lg->SetTextSize(0.06);
                    lg->Draw();
                }
                delete h; delete rp; delete hA; delete hB; delete h5;
            }
        bla::SaveFig(c, Form("%s/wp0_sim_mcerr", kFig));
        delete c;
    }
    fclose(o);
    f->Close();
    printf("wp0_check: done\n");
}
