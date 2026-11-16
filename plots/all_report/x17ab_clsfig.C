// ============================================================================
// x17ab_clsfig.C  --  WP9 X-b: figure and table of asymptotic against toy-based CLs (reads x17ab_toycls.C outputs)
// ============================================================================
// For every data set: toy CLs with binomial errors and the asymptotic CLs at the same r; the toy limit is the CLs = 0.10
// solution of a weighted quadratic fit of ln CLs in r (x17m::ToyUL). Table docs/latex/tables_all/x17ab_toycls.txt, figure
// docs/latex/figs_all/x17ab_toycls.pdf/.png (four panels: data sets 0, 3, 5, 7).
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_x17ab", kTRUE)' 'plots/all_report/x17ab_clsfig.C+O'
// ============================================================================
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"

#include "x17_model.h"

void x17ab_clsfig() {
    gStyle->SetOptStat(0);
    const char* dir = "/hdd3/przygoda/HADES/FAT_out/all_report/x17ab";
    FILE* tb = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/x17ab_toycls.txt", "w");
    fprintf(tb, "# x17ab_toycls.txt -- asymptotic against toy-based 90 %% CLs limits (x17ab_toycls.C, x17ab_clsfig.C); 1000 toys per r value and hypothesis;\n");
    fprintf(tb, "# toy limit from a weighted quadratic fit of ln CLs(r). Limits in BR_X BR_ee = r x 2.68e-4. Criterion of the spec: |toy/asymptotic - 1| <= 0.10.\n");
    fprintf(tb, "%-3s %-78s %11s %11s %7s\n", "set", "data set", "UL asympt.", "UL toys", "ratio");
    TCanvas c("c", "", 1400, 1100);
    c.Divide(2, 2);
    const int show[4] = {0, 3, 6, 7};
    int panel = 0;
    for (int w = 0; w <= 7; ++w) {
        std::ifstream in(Form("%s/toycls_%d.txt", dir, w));
        if (!in) continue;
        std::string line, label;
        std::vector<x17m::ToyCLsPoint> pts;
        std::vector<double> asy;
        double ulA = 0;
        while (std::getline(in, line)) {
            if (line.rfind("# x17ab_toycls.C data set", 0) == 0) { const size_t p = line.find(':'); label = line.substr(p + 2, line.find(';') - p - 2); continue; }
            if (line.empty() || line[0] == '#') continue;
            if (line.rfind("SUMMARY", 0) == 0) { std::istringstream ss(line); std::string a, b, d; ss >> a >> b >> d >> ulA; continue; }
            std::istringstream ss(line);
            x17m::ToyCLsPoint p; double br, bar_dummy; std::string bar;
            ss >> p.mu >> br >> p.qobs >> p.clsb >> p.clb >> p.cls >> p.clsErr >> bar >> bar_dummy;
            pts.push_back(p); asy.push_back(bar_dummy);
        }
        bool undefined = true;
        for (const auto& p : pts) if (p.clb > 0) undefined = false;
        const double ulT = undefined ? std::numeric_limits<double>::quiet_NaN() : x17m::ToyUL(pts);
        if (undefined)
            fprintf(tb, "%-3d %-78s %11.4e %11s %7s  (toy CLs undefined: q_obs above every toy under r and under 0, CL_b = 0)\n", w, label.c_str(), ulA * x17m::kBRdir, "-", "-");
        else
            fprintf(tb, "%-3d %-78s %11.4e %11.4e %7.3f\n", w, label.c_str(), ulA * x17m::kBRdir, ulT * x17m::kBRdir, ulT / ulA);
        printf("x17clsfig: set %d %s: asymptotic %.4e, toys %.4e, ratio %.3f\n", w, label.c_str(), ulA * x17m::kBRdir, ulT * x17m::kBRdir, ulT / ulA);
        bool drawIt = false;
        for (int q = 0; q < 4; ++q) if (show[q] == w) drawIt = true;
        if (!drawIt) continue;
        c.cd(++panel);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetLogy();
        TGraphErrors* gt = new TGraphErrors();
        TGraph* ga = new TGraph();
        for (size_t k = 0; k < pts.size(); ++k) {
            gt->SetPoint(k, pts[k].mu * x17m::kBRdir * 1e6, std::max(pts[k].cls, 1e-4)); gt->SetPointError(k, 0, pts[k].clsErr);
            ga->SetPoint(k, pts[k].mu * x17m::kBRdir * 1e6, asy[k]);
        }
        gt->SetTitle(Form("%s;r #times 2.68#times10^{-4} [10^{-6}];CL_{s}", label.c_str()));
        gt->SetMarkerStyle(20); gt->SetMinimum(3e-3); gt->SetMaximum(1.);
        gt->GetXaxis()->SetTitleSize(0.05); gt->GetYaxis()->SetTitleSize(0.05);
        gt->Draw("AP");
        ga->SetLineColor(kRed + 1); ga->SetLineWidth(2); ga->Draw("L SAME");
        const double xlo = gt->GetXaxis()->GetXmin(), xhi = gt->GetXaxis()->GetXmax();
        TLine* l = new TLine(xlo, 0.1, xhi, 0.1); l->SetLineStyle(2); l->Draw();
        TLegend* lg = new TLegend(0.17, 0.17, 0.62, 0.36); lg->SetBorderSize(0); lg->SetTextSize(0.04);
        lg->AddEntry(gt, Form("toys: UL %.2f", ulT * x17m::kBRdir * 1e6), "lp");
        lg->AddEntry(ga, Form("asymptotic: UL %.2f", ulA * x17m::kBRdir * 1e6), "l");
        lg->Draw();
    }
    fclose(tb);
    c.SaveAs("/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/x17ab_toycls.pdf");
    c.SaveAs("/home/przygoda/HADES/FAT/FAT/docs/latex/figs_all/x17ab_toycls.png");
}
