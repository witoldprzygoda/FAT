// research/plots/ecal_combined_map_view.C — visual inspection of the 3D
// ECAL energy-scale correction map produced by ecal_combined_3dscan.C.
//
// Reads:
//   ecal_combined_3dscan_<dataset>.root  with TH3D h_s_combined_3d
//   axes (E_γ, θ_γ, φ_γ_full)
//
// Produces (under plots/output/, suffix per dataset):
//   1) ecal_combined_map_view_<ds>_1d_E.{pdf,png}        — s vs E_γ
//                                                          (median over θ, φ)
//   2) ecal_combined_map_view_<ds>_1d_theta.{pdf,png}    — s vs θ_γ
//                                                          (median over E, φ)
//   3) ecal_combined_map_view_<ds>_1d_phi.{pdf,png}      — s vs φ_γ
//                                                          (median over E, θ)
//   4) ecal_combined_map_view_<ds>_2d_Eth.{pdf,png}      — s in (E, θ)
//                                                          (mean over φ)
//   5) ecal_combined_map_view_<ds>_2d_Ephi.{pdf,png}     — s in (E, φ)
//                                                          (mean over θ)
//   6) ecal_combined_map_view_<ds>_2d_thphi.{pdf,png}    — s in (θ, φ)
//                                                          (mean over E)
//   7) ecal_combined_map_view_<ds>_2d_thphi_perE.{pdf,png} — 2D (θ, φ)
//                                                          panels, one per
//                                                          E-bin (full map
//                                                          without averaging)
//
// Cells with content == 0 (empty or dead-sector) are excluded from medians/
// means; in 2D heatmaps they appear as white (zero suppressed via SetMinimum
// just below typical values, see kZMin/kZMax below).
//
// Usage (from research/):
//   root -l -b -q plots/ecal_combined_map_view.C            // exp default
//   root -l -b -q 'plots/ecal_combined_map_view.C("sim")'   // sim
//   root -l -b -q 'plots/ecal_combined_map_view.C("exp_oa4")' // any custom
//                                                           // suffix → reads
//                                                           // ecal_combined_3dscan_exp_oa4.root

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <utility>

namespace {
    // Color z-range for 2D heatmaps. The s map clusters around 1.0 with
    // few-percent variation; tightening the z-range makes structure visible.
    constexpr double kZMin = 0.85;
    constexpr double kZMax = 1.15;
    // Empty cells (zero content) → drawn as white (under min via kZMin).

    double medianOf(std::vector<double>& v) {
        if (v.empty()) return 0.0;
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
    }

    // Standard error of the median from a (small) sample, estimated as
    //   1.2533 * std(sample) / sqrt(N)
    // (the 1.2533 factor is the asymptotic ratio of σ_median/σ_mean for a
    // Gaussian; for non-Gaussian distributions this is a rough but
    // serviceable estimate of "how well the central value is determined").
    double medianStatErr(const std::vector<double>& v) {
        if (v.size() < 2) return 0.0;
        double s = 0, s2 = 0;
        for (double x : v) { s += x; s2 += x * x; }
        const double mean = s / v.size();
        const double var  = std::max(0.0, s2 / v.size() - mean * mean);
        const double std  = std::sqrt(var);
        return 1.2533 * std / std::sqrt((double)v.size());
    }

    void prettyAxes(TH1* h) {
        if (!h) return;
        h->GetXaxis()->SetTitleSize(0.045);
        h->GetYaxis()->SetTitleSize(0.045);
        h->GetXaxis()->SetLabelSize(0.040);
        h->GetYaxis()->SetLabelSize(0.040);
        h->GetXaxis()->SetTitleOffset(1.05);
        h->GetYaxis()->SetTitleOffset(1.05);
    }

    void drawSectorLines2DPhi(double y_lo, double y_hi) {
        // Vertical grid at sector boundaries (φ multiples of 60°). Used on
        // any 2D heatmap whose Y-axis is φ. Caller arranges that the canvas
        // is in the right state before drawing.
        for (int s = 1; s < 6; ++s) {
            auto* l = new TLine(y_lo, s * 60.0, y_hi, s * 60.0);
            l->SetLineStyle(3);
            l->SetLineColor(kGray + 2);
            l->Draw();
        }
    }

    void saveCanvas(TCanvas* c, const std::string& base) {
        c->SaveAs((base + ".pdf").c_str());
        c->SaveAs((base + ".png").c_str());
        std::cout << "  wrote " << base << ".{pdf,png}\n";
    }
}

void ecal_combined_map_view(const char* dataset = "exp") {

    const std::string ds   = dataset;
    const std::string fpath = "ecal_combined_3dscan_" + ds + ".root";
    TFile* f = TFile::Open(fpath.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << fpath << "\n"; return;
    }
    auto* h3 = (TH3D*)f->Get("h_s_combined_3d");
    if (!h3) {
        std::cerr << "TH3D 'h_s_combined_3d' not in " << fpath << "\n";
        f->Close(); return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gSystem->mkdir("plots/output", kTRUE);

    const int nE  = h3->GetNbinsX();
    const int nTh = h3->GetNbinsY();
    const int nPh = h3->GetNbinsZ();

    auto* axE  = h3->GetXaxis();
    auto* axTh = h3->GetYaxis();
    auto* axPh = h3->GetZaxis();

    // --------------------------------------------------------------------
    // Helper for 1D projection plots:
    // For each bin of the chosen axis, collect non-zero cell values along
    // the other two axes, compute median (central value) and a stat-error
    // estimate from the sample dispersion. Build a TGraphErrors with X
    // half-width as horizontal bar and the median stat-error as vertical
    // bar.
    // --------------------------------------------------------------------
    auto build1DGraph = [&](
        std::function<std::vector<double>(int)> collect,
        int n, std::function<double(int)> binCenter,
        std::function<double(int)> binHalfWidth) -> TGraphErrors*
    {
        std::vector<double> xs, exs, ys, eys;
        for (int i = 1; i <= n; ++i) {
            auto vs = collect(i);
            if (vs.empty()) continue;
            xs .push_back(binCenter(i));
            exs.push_back(binHalfWidth(i));
            ys .push_back(medianOf(vs));
            eys.push_back(medianStatErr(vs));
        }
        return new TGraphErrors((int)xs.size(),
                                xs.data(), ys.data(),
                                exs.data(), eys.data());
    };

    auto autoYRange = [&](TGraphErrors* g) {
        double y_lo = 1e9, y_hi = -1e9;
        for (int i = 0; i < g->GetN(); ++i) {
            const double y  = g->GetPointY(i);
            const double ey = g->GetErrorY(i);
            if (y - ey < y_lo) y_lo = y - ey;
            if (y + ey > y_hi) y_hi = y + ey;
        }
        // Pad ±20% of the data range, never tighter than ±2% around 1.
        const double dy = std::max(0.02, 0.2 * (y_hi - y_lo));
        return std::make_pair(y_lo - dy, y_hi + dy);
    };

    // --------------------------------------------------------------------
    // 1) 1D — s vs E (median over θ, φ), error bars + bin width
    // --------------------------------------------------------------------
    {
        auto* g = build1DGraph(
            [&](int e) {
                std::vector<double> vs;
                for (int t = 1; t <= nTh; ++t)
                  for (int p = 1; p <= nPh; ++p) {
                      const double v = h3->GetBinContent(e, t, p);
                      if (v > 0) vs.push_back(v);
                  }
                return vs;
            },
            nE,
            [&](int i){ return 0.5 * (axE->GetBinLowEdge(i) + axE->GetBinUpEdge(i)); },
            [&](int i){ return 0.5 * (axE->GetBinUpEdge(i)  - axE->GetBinLowEdge(i)); });
        g->SetTitle(";E_{#gamma} [GeV];median s(E)");
        g->SetMarkerStyle(20); g->SetMarkerSize(1.0); g->SetMarkerColor(kBlue+1);
        g->SetLineColor(kBlue+1); g->SetLineWidth(1);
        const auto yr = autoYRange(g);

        TCanvas c(("c1d_E_" + ds).c_str(), "s vs E", 1100, 700);
        c.SetMargin(0.13, 0.05, 0.14, 0.08); c.SetGrid();
        g->Draw("APE");
        g->GetXaxis()->SetLimits(axE->GetXmin(), axE->GetXmax());
        g->GetYaxis()->SetRangeUser(yr.first, yr.second);
        prettyAxes(g->GetHistogram());
        auto* l1 = new TLine(axE->GetXmin(), 1.0, axE->GetXmax(), 1.0);
        l1->SetLineStyle(2); l1->SetLineColor(kRed); l1->SetLineWidth(2); l1->Draw();
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
        lat.DrawLatex(0.15, 0.92,
            TString::Format("dataset = %s   (median over #theta, #phi cells; X-bar = bin width, Y-bar = stat. err on median)",
                            ds.c_str()).Data());
        saveCanvas(&c, "plots/output/ecal_combined_map_view_" + ds + "_1d_E");
    }

    // --------------------------------------------------------------------
    // 2) 1D — s vs θ (median over E, φ)
    // --------------------------------------------------------------------
    {
        auto* g = build1DGraph(
            [&](int t) {
                std::vector<double> vs;
                for (int e = 1; e <= nE; ++e)
                  for (int p = 1; p <= nPh; ++p) {
                      const double v = h3->GetBinContent(e, t, p);
                      if (v > 0) vs.push_back(v);
                  }
                return vs;
            },
            nTh,
            [&](int i){ return 0.5 * (axTh->GetBinLowEdge(i) + axTh->GetBinUpEdge(i)); },
            [&](int i){ return 0.5 * (axTh->GetBinUpEdge(i)  - axTh->GetBinLowEdge(i)); });
        g->SetTitle(";#theta_{#gamma} [deg];median s(#theta)");
        g->SetMarkerStyle(20); g->SetMarkerSize(1.0); g->SetMarkerColor(kBlue+1);
        g->SetLineColor(kBlue+1); g->SetLineWidth(1);
        const auto yr = autoYRange(g);

        TCanvas c(("c1d_theta_" + ds).c_str(), "s vs theta", 1100, 700);
        c.SetMargin(0.13, 0.05, 0.14, 0.08); c.SetGrid();
        g->Draw("APE");
        g->GetXaxis()->SetLimits(axTh->GetXmin(), axTh->GetXmax());
        g->GetYaxis()->SetRangeUser(yr.first, yr.second);
        prettyAxes(g->GetHistogram());
        auto* l1 = new TLine(axTh->GetXmin(), 1.0, axTh->GetXmax(), 1.0);
        l1->SetLineStyle(2); l1->SetLineColor(kRed); l1->SetLineWidth(2); l1->Draw();
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
        lat.DrawLatex(0.15, 0.92,
            TString::Format("dataset = %s   (median over E_{#gamma}, #phi cells; X-bar = bin width, Y-bar = stat. err on median)",
                            ds.c_str()).Data());
        saveCanvas(&c, "plots/output/ecal_combined_map_view_" + ds + "_1d_theta");
    }

    // --------------------------------------------------------------------
    // 3) 1D — s vs φ_full (median over E, θ)
    // --------------------------------------------------------------------
    {
        auto* g = build1DGraph(
            [&](int p) {
                std::vector<double> vs;
                for (int e = 1; e <= nE; ++e)
                  for (int t = 1; t <= nTh; ++t) {
                      const double v = h3->GetBinContent(e, t, p);
                      if (v > 0) vs.push_back(v);
                  }
                return vs;
            },
            nPh,
            [&](int i){ return 0.5 * (axPh->GetBinLowEdge(i) + axPh->GetBinUpEdge(i)); },
            [&](int i){ return 0.5 * (axPh->GetBinUpEdge(i)  - axPh->GetBinLowEdge(i)); });
        g->SetTitle(";#phi_{#gamma} [deg];median s(#phi)");
        g->SetMarkerStyle(20); g->SetMarkerSize(0.9); g->SetMarkerColor(kBlack);
        g->SetLineColor(kBlack); g->SetLineWidth(1);
        const auto yr = autoYRange(g);

        TCanvas c(("c1d_phi_" + ds).c_str(), "s vs phi", 1300, 600);
        c.SetMargin(0.10, 0.05, 0.13, 0.08); c.SetGrid();
        g->Draw("APE");
        g->GetXaxis()->SetLimits(0, 360);
        g->GetYaxis()->SetRangeUser(yr.first, yr.second);
        prettyAxes(g->GetHistogram());
        for (int s = 1; s < 6; ++s) {
            auto* l = new TLine(s * 60.0, yr.first, s * 60.0, yr.second);
            l->SetLineStyle(3); l->SetLineColor(kGray + 2); l->Draw();
        }
        auto* l1 = new TLine(0, 1.0, 360, 1.0);
        l1->SetLineStyle(2); l1->SetLineColor(kRed); l1->SetLineWidth(2); l1->Draw();
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
        lat.DrawLatex(0.12, 0.92,
            TString::Format("dataset = %s   (median over E_{#gamma}, #theta cells; X-bar = bin width, Y-bar = stat. err on median)",
                            ds.c_str()).Data());
        saveCanvas(&c, "plots/output/ecal_combined_map_view_" + ds + "_1d_phi");
    }

    // --------------------------------------------------------------------
    // 4) 2D — s in (E, θ), mean over φ (excluding empty cells)
    // --------------------------------------------------------------------
    {
        TH2D* h = new TH2D(("h_s_Eth_" + ds).c_str(),
            ";E_{#gamma} [GeV];#theta_{#gamma} [deg];mean s   (over #phi)",
            nE,  axE ->GetXbins()->GetArray(),
            nTh, axTh->GetXbins()->GetArray());
        for (int e = 1; e <= nE; ++e)
          for (int t = 1; t <= nTh; ++t) {
              double sum = 0; int n = 0;
              for (int p = 1; p <= nPh; ++p) {
                  const double v = h3->GetBinContent(e, t, p);
                  if (v > 0) { sum += v; ++n; }
              }
              h->SetBinContent(e, t, n > 0 ? sum / n : 0);
          }
        h->SetMinimum(kZMin); h->SetMaximum(kZMax);
        prettyAxes(h);
        TCanvas c(("c2d_Eth_" + ds).c_str(), "s(E,theta)", 1100, 700);
        c.SetMargin(0.12, 0.14, 0.13, 0.08); c.SetGrid();
        h->Draw("COLZ");
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.04);
        lat.DrawLatex(0.13, 0.92,
            TString::Format("dataset = %s", ds.c_str()).Data());
        saveCanvas(&c, "plots/output/ecal_combined_map_view_" + ds + "_2d_Eth");
    }

    // --------------------------------------------------------------------
    // 5) 2D — s in (E, φ), mean over θ
    // --------------------------------------------------------------------
    {
        TH2D* h = new TH2D(("h_s_Ephi_" + ds).c_str(),
            ";E_{#gamma} [GeV];#phi_{#gamma} [deg];mean s   (over #theta)",
            nE, axE->GetXbins()->GetArray(),
            nPh, axPh->GetXmin(), axPh->GetXmax());
        for (int e = 1; e <= nE; ++e)
          for (int p = 1; p <= nPh; ++p) {
              double sum = 0; int n = 0;
              for (int t = 1; t <= nTh; ++t) {
                  const double v = h3->GetBinContent(e, t, p);
                  if (v > 0) { sum += v; ++n; }
              }
              h->SetBinContent(e, p, n > 0 ? sum / n : 0);
          }
        h->SetMinimum(kZMin); h->SetMaximum(kZMax);
        prettyAxes(h);
        TCanvas c(("c2d_Ephi_" + ds).c_str(), "s(E,phi)", 1100, 700);
        c.SetMargin(0.12, 0.14, 0.13, 0.08); c.SetGrid();
        h->Draw("COLZ");
        // Sector boundaries on phi axis (Y).
        drawSectorLines2DPhi(axE->GetXmin(), axE->GetXmax());
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.04);
        lat.DrawLatex(0.13, 0.92,
            TString::Format("dataset = %s", ds.c_str()).Data());
        saveCanvas(&c, "plots/output/ecal_combined_map_view_" + ds + "_2d_Ephi");
    }

    // --------------------------------------------------------------------
    // 6) 2D — s in (θ, φ), mean over E
    // --------------------------------------------------------------------
    {
        TH2D* h = new TH2D(("h_s_thphi_" + ds).c_str(),
            ";#theta_{#gamma} [deg];#phi_{#gamma} [deg];mean s   (over E_{#gamma})",
            nTh, axTh->GetXbins()->GetArray(),
            nPh, axPh->GetXmin(), axPh->GetXmax());
        for (int t = 1; t <= nTh; ++t)
          for (int p = 1; p <= nPh; ++p) {
              double sum = 0; int n = 0;
              for (int e = 1; e <= nE; ++e) {
                  const double v = h3->GetBinContent(e, t, p);
                  if (v > 0) { sum += v; ++n; }
              }
              h->SetBinContent(t, p, n > 0 ? sum / n : 0);
          }
        h->SetMinimum(kZMin); h->SetMaximum(kZMax);
        prettyAxes(h);
        TCanvas c(("c2d_thphi_" + ds).c_str(), "s(theta,phi)", 1100, 700);
        c.SetMargin(0.12, 0.14, 0.13, 0.08); c.SetGrid();
        h->Draw("COLZ");
        drawSectorLines2DPhi(axTh->GetXmin(), axTh->GetXmax());
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.04);
        lat.DrawLatex(0.13, 0.92,
            TString::Format("dataset = %s", ds.c_str()).Data());
        saveCanvas(&c, "plots/output/ecal_combined_map_view_" + ds + "_2d_thphi");
    }

    // --------------------------------------------------------------------
    // 7) 2D panels (θ, φ), one per E-bin — full map without averaging
    // --------------------------------------------------------------------
    {
        const int nx = (nE <= 2) ? nE : 2;
        const int ny = (nE + nx - 1) / nx;
        TCanvas c(("c2d_thphi_perE_" + ds).c_str(),
                  "s(theta,phi) per E-bin",
                  600 * nx, 450 * ny);
        c.Divide(nx, ny);
        for (int e = 1; e <= nE; ++e) {
            c.cd(e);
            gPad->SetMargin(0.12, 0.14, 0.13, 0.10);
            gPad->SetGrid();
            const TString hname =
                TString::Format("h_s_thphi_E%d_%s", e, ds.c_str());
            TH2D* h = new TH2D(hname,
                TString::Format(
                    "E_{#gamma} #in [%.2f, %.2f] GeV;"
                    "#theta_{#gamma} [deg];#phi_{#gamma} [deg]",
                    axE->GetBinLowEdge(e), axE->GetBinUpEdge(e)),
                nTh, axTh->GetXbins()->GetArray(),
                nPh, axPh->GetXmin(), axPh->GetXmax());
            for (int t = 1; t <= nTh; ++t)
              for (int p = 1; p <= nPh; ++p) {
                  const double v = h3->GetBinContent(e, t, p);
                  h->SetBinContent(t, p, v);  // 0 → empty (below kZMin)
              }
            h->SetMinimum(kZMin); h->SetMaximum(kZMax);
            prettyAxes(h);
            h->Draw("COLZ");
            drawSectorLines2DPhi(axTh->GetXmin(), axTh->GetXmax());
        }
        saveCanvas(&c, "plots/output/ecal_combined_map_view_" + ds + "_2d_thphi_perE");
    }

    // --------------------------------------------------------------------
    // 8) 1D — s vs E, **per φ bin** (slices). 36 panels in 6×6 grid:
    //    rows = 6 sectors, columns = 6 φ-within-sector bins.
    //    Each panel shows s(E) for one fixed φ slice, median over θ cells.
    //    Useful to assess whether s(E) per φ has a smooth trend that
    //    could be parametrized analytically.
    // --------------------------------------------------------------------
    {
        TCanvas c(("c1d_E_perphi_" + ds).c_str(),
                  "s vs E per phi slice", 1900, 1500);
        c.Divide(6, 6, 0.001, 0.001);
        for (int p = 1; p <= nPh; ++p) {
            c.cd(p);
            gPad->SetMargin(0.18, 0.05, 0.20, 0.10);
            gPad->SetGrid();

            std::vector<double> xs, exs, ys, eys;
            for (int e = 1; e <= nE; ++e) {
                std::vector<double> vs;
                for (int t = 1; t <= nTh; ++t) {
                    const double v = h3->GetBinContent(e, t, p);
                    if (v > 0) vs.push_back(v);
                }
                if (vs.empty()) continue;
                xs .push_back(0.5 * (axE->GetBinLowEdge(e) + axE->GetBinUpEdge(e)));
                exs.push_back(0.5 * (axE->GetBinUpEdge(e)  - axE->GetBinLowEdge(e)));
                ys .push_back(medianOf(vs));
                eys.push_back(medianStatErr(vs));
            }
            if (xs.empty()) {
                // Empty (e.g. dead sector 1) — draw stub.
                TLatex lat; lat.SetNDC(); lat.SetTextSize(0.13);
                lat.DrawLatex(0.20, 0.45, TString::Format(
                    "#phi #in [%.0f, %.0f) (empty)",
                    axPh->GetBinLowEdge(p), axPh->GetBinUpEdge(p)).Data());
                continue;
            }
            auto* g = new TGraphErrors((int)xs.size(),
                xs.data(), ys.data(), exs.data(), eys.data());
            g->SetTitle(TString::Format(
                "#phi #in [%.0f, %.0f);E_{#gamma} [GeV];s",
                axPh->GetBinLowEdge(p), axPh->GetBinUpEdge(p)));
            g->SetMarkerStyle(20); g->SetMarkerSize(0.7);
            const int sec = (p - 1) / 6;  // 0..5
            const int col = (sec == 0) ? kBlue+1
                          : (sec == 1) ? kRed+1
                          : (sec == 2) ? kGreen+2
                          : (sec == 3) ? kMagenta+1
                          : (sec == 4) ? kOrange+7
                                       : kCyan+2;
            g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(1);
            g->Draw("APE");
            g->GetXaxis()->SetLimits(axE->GetXmin(), axE->GetXmax());
            g->GetYaxis()->SetRangeUser(kZMin, kZMax);
            g->GetXaxis()->SetTitleSize(0.08);
            g->GetYaxis()->SetTitleSize(0.08);
            g->GetXaxis()->SetLabelSize(0.07);
            g->GetYaxis()->SetLabelSize(0.07);
            g->GetXaxis()->SetTitleOffset(1.1);
            g->GetYaxis()->SetTitleOffset(0.9);
            // PDG line at s=1 for reference.
            auto* l1 = new TLine(axE->GetXmin(), 1.0, axE->GetXmax(), 1.0);
            l1->SetLineStyle(2); l1->SetLineColor(kBlack); l1->SetLineWidth(1);
            l1->Draw();
        }
        saveCanvas(&c, "plots/output/ecal_combined_map_view_" + ds + "_1d_E_per_phi");
    }

    // --------------------------------------------------------------------
    // 9) 1D — s vs E, **per θ bin** (slices). nθ panels in (4 × ceil(nθ/4)).
    //    Each panel shows s(E) for one fixed θ slice, median over φ cells.
    // --------------------------------------------------------------------
    {
        const int nx = 4;
        const int ny = (nTh + nx - 1) / nx;
        TCanvas c(("c1d_E_pertheta_" + ds).c_str(),
                  "s vs E per theta slice", 600 * nx, 400 * ny);
        c.Divide(nx, ny);
        for (int t = 1; t <= nTh; ++t) {
            c.cd(t);
            gPad->SetMargin(0.13, 0.05, 0.14, 0.10);
            gPad->SetGrid();

            std::vector<double> xs, exs, ys, eys;
            for (int e = 1; e <= nE; ++e) {
                std::vector<double> vs;
                for (int p = 1; p <= nPh; ++p) {
                    const double v = h3->GetBinContent(e, t, p);
                    if (v > 0) vs.push_back(v);
                }
                if (vs.empty()) continue;
                xs .push_back(0.5 * (axE->GetBinLowEdge(e) + axE->GetBinUpEdge(e)));
                exs.push_back(0.5 * (axE->GetBinUpEdge(e)  - axE->GetBinLowEdge(e)));
                ys .push_back(medianOf(vs));
                eys.push_back(medianStatErr(vs));
            }
            if (xs.empty()) continue;
            auto* g = new TGraphErrors((int)xs.size(),
                xs.data(), ys.data(), exs.data(), eys.data());
            g->SetTitle(TString::Format(
                "#theta #in [%.0f, %.0f) deg;E_{#gamma} [GeV];median s",
                axTh->GetBinLowEdge(t), axTh->GetBinUpEdge(t)));
            g->SetMarkerStyle(20); g->SetMarkerSize(1.0); g->SetMarkerColor(kBlue+1);
            g->SetLineColor(kBlue+1); g->SetLineWidth(1);
            g->Draw("APE");
            g->GetXaxis()->SetLimits(axE->GetXmin(), axE->GetXmax());
            g->GetYaxis()->SetRangeUser(kZMin, kZMax);
            auto* l1 = new TLine(axE->GetXmin(), 1.0, axE->GetXmax(), 1.0);
            l1->SetLineStyle(2); l1->SetLineColor(kRed); l1->SetLineWidth(2);
            l1->Draw();
        }
        saveCanvas(&c, "plots/output/ecal_combined_map_view_" + ds + "_1d_E_per_theta");
    }

    f->Close();
    std::cout << "\nDone. 9 plots in plots/output/ecal_combined_map_view_"
              << ds << "_*.{pdf,png}\n";
}
