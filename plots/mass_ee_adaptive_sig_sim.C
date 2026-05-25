// mass_ee_adaptive_sig_sim.C — adaptive-binning M_ee for SMASH simulation
// with trigger-flag splitting. Analog of mass_ee_adaptive_sig_exp.C but:
//   * no CB subtraction (sim is pure signal after MC purity gate)
//   * no factor 63 (sim has no PT2 downscale — every event evaluated for both
//     triggers independently)
//   * three histogram variants: H (no trigger), H_pt3, H_pt2
//
// Two canvases:
//   1) <out>_adaptive_sig.{pdf,png}            — three panels:
//        • spectra overlay (none / PT3 / PT2)               log Y
//        • trigger efficiencies (PT3/none, PT2/none)        linear Y, 0–1.2
//        • fine ratio N_PT2 / N_PT3                          linear Y
//   2) <out>_adaptive_sig_ratios_coarse.{pdf,png} — three coarser panels:
//        • N_PT2 / N_PT3                       (correction factor)
//        • N_PT3 / N_PT2                       (PT3-conditional-on-PT2 eff.)
//        • N_PT3 / N_none                      (absolute PT3 efficiency)
//
// Adaptive binning: PT3-statistics-driven (PT3 is the most-filtered sample
// so it sets the floor). Walk a fine 1 MeV grid; close each bin when the
// accumulated PT3 weight reaches `N_min_pt3_spectra` OR width hits
// `max_bin_width`. The coarser ratio canvas uses `N_min_pt3_ratio`.
//
// All fills are weighted by `sim_genweight` (the per-event SMASH luminosity
// weight stored on every dilepton_nt row).
//
// Usage:
//   root -l -b -q plots/mass_ee_adaptive_sig_sim.C
//   root -l -b -q 'plots/mass_ee_adaptive_sig_sim.C(500)'              # wider spectra
//   root -l -b -q 'plots/mass_ee_adaptive_sig_sim.C(200, 0.2, 5000)'    # coarser ratios

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr const char* kInputFile = "output_epem_sim.root";
constexpr const char* kNtName    = "dilepton_nt";
constexpr const char* kVarExpr   = "m_ee";

// Detect whether sim_genweight is filled in the loaded ntuple. SMASH
// lepton mode stores the dilepton-process weight there; SMASH std mode
// leaves it at 0. When the global integral is zero we fall back to
// unweighted ("1.0"); otherwise we use the per-event weight. Returns a
// TTree::Draw weight expression as a std::string. Macro-level decision
// — NOT per-event — so a few stray sim_genweight==0 events in lepton
// data don't suddenly get over-weighted by a fallback ternary.
std::string detectWeightExpr(TTree* nt) {
    TH1D* h = new TH1D("h_wsum_detect", "", 1, 0, 2);
    nt->Draw("1>>h_wsum_detect", "sim_genweight", "goff");
    const double sum = h->Integral();
    delete h;
    if (sum > 0.0) {
        std::cout << "[sim] sim_genweight integral = " << sum
                  << "  → using weighted fills (sim_genweight)\n";
        return "sim_genweight";
    } else {
        std::cout << "[sim] sim_genweight integral = 0  → unweighted (1.0)\n";
        return "1.0";
    }
}
constexpr double      kXmin      = 0.0;
constexpr double      kXmax      = 1.4;
constexpr int         kFineNb    = 1400;       // 1 MeV/bin seed

TH1D* drawFromNt(TTree* t, const std::string& cut_w, const std::string& name,
                 int nb, const double* edges) {
    TH1D* h = new TH1D(name.c_str(), "", nb, edges);
    h->Sumw2();
    t->Draw((std::string(kVarExpr) + ">>" + name).c_str(), cut_w.c_str(), "goff");
    h->SetDirectory(nullptr);
    return h;
}

TH1D* makeRatio(TH1D* num, TH1D* den, const std::string& name,
                double scale = 1.0) {
    TH1D* r = static_cast<TH1D*>(num->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    const int nb = num->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double n  = num->GetBinContent(b);
        const double en = num->GetBinError(b);
        const double d  = den->GetBinContent(b);
        const double ed = den->GetBinError(b);
        if (d != 0.0 && std::isfinite(d) && std::isfinite(n)) {
            const double val = scale * n / d;
            const double rel = std::sqrt(
                (n != 0.0 ? std::pow(en / n, 2) : 0.0) +
                std::pow(ed / d, 2));
            r->SetBinContent(b, val);
            r->SetBinError  (b, std::abs(val) * rel);
        }
    }
    return r;
}

double safeMaxAbs(TH1D* h) {
    double mx = 0.0;
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = std::abs(h->GetBinContent(b)) + h->GetBinError(b);
        if (std::isfinite(v) && v > mx) mx = v;
    }
    return mx;
}

double safeMaxRobust(TH1D* h, double rel_err_cap = 0.50) {
    double mx = 0.0;
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = h->GetBinContent(b);
        const double e = h->GetBinError(b);
        if (!std::isfinite(v) || v <= 0.0) continue;
        if (e > 0.0 && e / std::abs(v) > rel_err_cap) continue;
        if (v + e > mx) mx = v + e;
    }
    return mx;
}

double safeMinPositive(TH1D* h) {
    double mn = std::numeric_limits<double>::infinity();
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double v = h->GetBinContent(b);
        if (std::isfinite(v) && v > 0.0 && v < mn) mn = v;
    }
    return std::isfinite(mn) ? mn : 1.0;
}

void styleNone(TH1D* h) { h->SetMarkerStyle(20); h->SetMarkerSize(0.65);
                          h->SetMarkerColor(kBlack); h->SetLineColor(kBlack); }
void stylePT3 (TH1D* h) { h->SetMarkerStyle(21); h->SetMarkerSize(0.65);
                          h->SetMarkerColor(kBlue + 1); h->SetLineColor(kBlue + 1); }
void stylePT2 (TH1D* h) { h->SetMarkerStyle(22); h->SetMarkerSize(0.7);
                          h->SetMarkerColor(kRed + 1); h->SetLineColor(kRed + 1); }

void styleAxes(TH1D* h) {
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

/// Adaptive edges from a weighted PT3 seed histogram. `w_expr` is the
/// per-event weight (`sim_genweight` for lepton sample, `1.0` for std).
std::vector<double> deriveEdges(TTree* nt, const std::string& w_expr,
                                double N_min, double max_bw,
                                const std::string& seed_name) {
    const double fine_w = (kXmax - kXmin) / kFineNb;
    TH1D* h_seed = new TH1D(seed_name.c_str(), "", kFineNb, kXmin, kXmax);
    h_seed->Sumw2();
    // PT3 sample drives the schedule (most-selective sample).
    nt->Draw((std::string(kVarExpr) + ">>" + seed_name).c_str(),
             (std::string("(pt3==1) * ") + w_expr).c_str(), "goff");
    h_seed->SetDirectory(nullptr);

    std::vector<double> edges;
    edges.reserve(kFineNb + 1);
    edges.push_back(kXmin);
    double acc = 0.0;
    for (int b = 1; b <= kFineNb; ++b) {
        acc += h_seed->GetBinContent(b);
        const double upper = h_seed->GetBinLowEdge(b + 1);
        const double width = upper - edges.back();
        if (acc >= N_min || width >= max_bw - 0.5 * fine_w) {
            edges.push_back(upper);
            acc = 0.0;
        }
    }
    if (edges.back() < kXmax - 0.5 * fine_w) edges.push_back(kXmax);
    delete h_seed;
    return edges;
}

}  // anonymous namespace

void mass_ee_adaptive_sig_sim(double N_min_pt3_spectra = 200.0,
                              double max_bin_width      = 0.2,
                              double N_min_pt3_ratio    = 2000.0) {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile* f = TFile::Open(kInputFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << kInputFile << "\n";
        return;
    }
    TTree* nt = dynamic_cast<TTree*>(f->Get(kNtName));
    if (!nt) {
        std::cerr << "Tree '" << kNtName << "' missing\n";
        return;
    }

    // Detect once whether weights are filled (lepton mode) or not (std mode).
    const std::string w_expr = detectWeightExpr(nt);

    // --------------------------------------------------------------------
    // FINE binning for spectra
    // --------------------------------------------------------------------
    const std::vector<double> edges = deriveEdges(
        nt, w_expr, N_min_pt3_spectra, max_bin_width, "h_seed_pt3_sim_fine");
    const int nb = static_cast<int>(edges.size()) - 1;
    std::cout << "[sim] spectra: " << nb << " bins, "
              << "N_min_pt3=" << N_min_pt3_spectra
              << ", max_width=" << max_bin_width << "\n";

    const std::string w_none = w_expr;
    const std::string w_pt3  = std::string("(pt3==1) * ") + w_expr;
    const std::string w_pt2  = std::string("(pt2==1) * ") + w_expr;

    TH1D* h_none = drawFromNt(nt, w_none, "h_none", nb, edges.data());
    TH1D* h_pt3  = drawFromNt(nt, w_pt3,  "h_pt3",  nb, edges.data());
    TH1D* h_pt2  = drawFromNt(nt, w_pt2,  "h_pt2",  nb, edges.data());

    // Fine ratios (signal-equivalent: in sim there's no CB to subtract).
    TH1D* r_pt2_pt3_fine = makeRatio(h_pt2, h_pt3, "r_pt2_pt3_fine");
    TH1D* eff_pt3_fine   = makeRatio(h_pt3, h_none, "eff_pt3_fine");
    TH1D* eff_pt2_fine   = makeRatio(h_pt2, h_none, "eff_pt2_fine");

    // Scale spectra by bin width — counts / GeV/c² (avoid spurious steps
    // from variable-width binning).
    for (TH1D* h : {h_none, h_pt3, h_pt2}) h->Scale(1.0, "width");

    styleNone(h_none); stylePT3(h_pt3); stylePT2(h_pt2);
    stylePT3(eff_pt3_fine); stylePT2(eff_pt2_fine);
    stylePT3(r_pt2_pt3_fine);
    for (TH1D* h : {h_none, h_pt3, h_pt2,
                    r_pt2_pt3_fine, eff_pt3_fine, eff_pt2_fine}) styleAxes(h);

    h_none->SetTitle("M_{ee} spectra (none / PT3 / PT2);"
                     "M_{e^{+}e^{-}} [GeV/c^{2}];"
                     "Weighted counts / (GeV/c^{2})");
    eff_pt3_fine->SetTitle("Trigger efficiency vs unbiased;"
                           "M_{e^{+}e^{-}} [GeV/c^{2}];"
                           "N_{trig} / N_{none}");
    r_pt2_pt3_fine->SetTitle("N_{PT2} / N_{PT3} (signal, fine bins);"
                             "M_{e^{+}e^{-}} [GeV/c^{2}];"
                             "N_{PT2} / N_{PT3}");

    // --------------------------------------------------------------------
    // Canvas 1: spectra overlay, efficiencies, fine ratio
    // --------------------------------------------------------------------
    TCanvas* c = new TCanvas("c_mass_ee_adaptive_sig_sim",
                             "SMASH M_{ee}: spectra / efficiencies / ratio",
                             1800, 600);
    c->Divide(3, 1, 0.001, 0.001);

    // Pad 1: spectra overlay (log Y)
    c->cd(1);
    gPad->SetLogy();
    gPad->SetMargin(0.14, 0.04, 0.13, 0.10);
    {
        const double ymax = std::max({safeMaxAbs(h_none), safeMaxAbs(h_pt3),
                                      safeMaxAbs(h_pt2)});
        const double ymin = std::max(0.5, std::min({safeMinPositive(h_none),
                                                    safeMinPositive(h_pt3),
                                                    safeMinPositive(h_pt2)}));
        h_none->GetYaxis()->SetRangeUser(ymin * 0.5, ymax * 4.0);
    }
    h_none->Draw("E1");
    h_pt3 ->Draw("E1 SAME");
    h_pt2 ->Draw("E1 SAME");
    {
        TLegend* leg = new TLegend(0.55, 0.74, 0.95, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
        leg->AddEntry(h_none, "no trigger (unbiased)", "lpe");
        leg->AddEntry(h_pt3,  "PT3 (bit 13)",          "lpe");
        leg->AddEntry(h_pt2,  "PT2 (bit 12)",          "lpe");
        leg->Draw();
    }

    // Pad 2: trigger efficiencies (linear Y, ~0–1)
    c->cd(2);
    gPad->SetMargin(0.14, 0.04, 0.13, 0.10);
    {
        const double ymax = std::max(safeMaxRobust(eff_pt3_fine),
                                     safeMaxRobust(eff_pt2_fine));
        const double y_hi = std::max(1.2, ymax * 1.10);
        eff_pt3_fine->SetTitle("Trigger efficiency vs unbiased;"
                               "M_{e^{+}e^{-}} [GeV/c^{2}];"
                               "N_{trig} / N_{none}");
        eff_pt3_fine->GetYaxis()->SetRangeUser(0.0, y_hi);
    }
    eff_pt3_fine->Draw("E1");
    eff_pt2_fine->Draw("E1 SAME");
    TLine* l_one = new TLine(kXmin, 1.0, kXmax, 1.0);
    l_one->SetLineStyle(3); l_one->SetLineColor(kGray + 2); l_one->Draw();
    {
        TLegend* leg = new TLegend(0.55, 0.74, 0.95, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.038);
        leg->AddEntry(eff_pt3_fine, "#epsilon_{PT3} = N_{PT3}/N_{none}", "lpe");
        leg->AddEntry(eff_pt2_fine, "#epsilon_{PT2} = N_{PT2}/N_{none}", "lpe");
        leg->AddEntry(l_one,        "1.0 (catches all)",                  "l");
        leg->Draw();
    }

    // Pad 3: fine ratio PT2/PT3
    c->cd(3);
    gPad->SetMargin(0.14, 0.04, 0.13, 0.10);
    {
        const double ymax = safeMaxRobust(r_pt2_pt3_fine);
        const double y_hi = (ymax > 0.0 && std::isfinite(ymax)) ? ymax * 1.25 : 2.0;
        r_pt2_pt3_fine->GetYaxis()->SetRangeUser(0.0, y_hi);
    }
    r_pt2_pt3_fine->Draw("E1");
    TLine* l_one_r = new TLine(kXmin, 1.0, kXmax, 1.0);
    l_one_r->SetLineStyle(3); l_one_r->SetLineColor(kGray + 2); l_one_r->Draw();
    {
        TLegend* leg = new TLegend(0.50, 0.78, 0.95, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.038);
        leg->AddEntry(r_pt2_pt3_fine, "N_{PT2} / N_{PT3}",   "lpe");
        leg->AddEntry(l_one_r,        "1.0 (PT3=PT2)",       "l");
        leg->Draw();
    }

    gSystem->mkdir("plots/output", true);
    c->SaveAs("plots/output/mass_ee_adaptive_sig_sim.pdf");
    c->SaveAs("plots/output/mass_ee_adaptive_sig_sim.png");
    std::cout << "Saved: plots/output/mass_ee_adaptive_sig_sim.{pdf,png}\n";

    // --------------------------------------------------------------------
    // Canvas 2: coarse ratios (signal only, three views)
    // --------------------------------------------------------------------
    const std::vector<double> edges_c = deriveEdges(
        nt, w_expr, N_min_pt3_ratio, max_bin_width, "h_seed_pt3_sim_coarse");
    const int nb_c = static_cast<int>(edges_c.size()) - 1;
    std::cout << "[sim] coarse: " << nb_c << " bins, "
              << "N_min_pt3_ratio=" << N_min_pt3_ratio << "\n";

    TH1D* h_none_c = drawFromNt(nt, w_none, "h_none_c", nb_c, edges_c.data());
    TH1D* h_pt3_c  = drawFromNt(nt, w_pt3,  "h_pt3_c",  nb_c, edges_c.data());
    TH1D* h_pt2_c  = drawFromNt(nt, w_pt2,  "h_pt2_c",  nb_c, edges_c.data());

    TH1D* r_pt2_pt3_c  = makeRatio(h_pt2_c, h_pt3_c, "r_pt2_pt3_c");
    TH1D* r_pt3_pt2_c  = makeRatio(h_pt3_c, h_pt2_c, "r_pt3_pt2_c");
    TH1D* r_pt3_none_c = makeRatio(h_pt3_c, h_none_c, "r_pt3_none_c");
    stylePT3(r_pt2_pt3_c); stylePT3(r_pt3_pt2_c); stylePT3(r_pt3_none_c);
    styleAxes(r_pt2_pt3_c); styleAxes(r_pt3_pt2_c); styleAxes(r_pt3_none_c);

    TCanvas* c2 = new TCanvas(
        "c_mass_ee_adaptive_sig_sim_ratios_coarse",
        "SMASH M_{ee} coarse ratios: PT2/PT3, PT3/PT2, PT3/none",
        2100, 600);
    c2->Divide(3, 1, 0.001, 0.001);

    auto drawRatioPad = [&](int idx, TH1D* h, const std::string& title,
                            const std::string& ylab, const std::string& leg_label,
                            double ref_y, const std::string& ref_label,
                            int ref_style, int ref_color, double y_max_floor) {
        c2->cd(idx);
        gPad->SetMargin(0.14, 0.04, 0.13, 0.10);
        const double ymax = safeMaxRobust(h);
        const double y_hi = (ymax > 0.0 && std::isfinite(ymax)) ? ymax * 1.25 : y_max_floor;
        h->SetTitle((title + ";M_{e^{+}e^{-}} [GeV/c^{2}];" + ylab).c_str());
        h->GetYaxis()->SetRangeUser(0.0, std::max(y_hi, y_max_floor));
        h->Draw("E1");
        TLine* lref = new TLine(kXmin, ref_y, kXmax, ref_y);
        lref->SetLineStyle(ref_style); lref->SetLineColor(ref_color);
        if (ref_style == 2) lref->SetLineWidth(2);
        lref->Draw();
        TLegend* leg = new TLegend(0.45, 0.78, 0.95, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.038);
        leg->AddEntry(h,    leg_label.c_str(), "lpe");
        leg->AddEntry(lref, ref_label.c_str(), "l");
        leg->Draw();
    };

    drawRatioPad(1, r_pt2_pt3_c,
        "N_{PT2}/N_{PT3} (coarse bins)",
        "N_{PT2}/N_{PT3}",
        "N_{PT2}/N_{PT3}  (= 1/#epsilon_{PT3|PT2})",
        1.0, "1.0 (PT3 = PT2)", 3, kGray + 2, 2.0);

    drawRatioPad(2, r_pt3_pt2_c,
        "N_{PT3}/N_{PT2} (coarse bins)",
        "N_{PT3}/N_{PT2}",
        "N_{PT3}/N_{PT2}  (= #epsilon_{PT3|PT2})",
        1.0, "1.0 (PT3 = PT2)", 3, kGray + 2, 1.2);

    drawRatioPad(3, r_pt3_none_c,
        "N_{PT3}/N_{none} (coarse bins)",
        "N_{PT3}/N_{none}",
        "N_{PT3}/N_{none}  (#epsilon_{PT3} vs unbiased)",
        1.0, "1.0 (PT3 catches all)", 3, kGray + 2, 1.2);

    c2->SaveAs("plots/output/mass_ee_adaptive_sig_sim_ratios_coarse.pdf");
    c2->SaveAs("plots/output/mass_ee_adaptive_sig_sim_ratios_coarse.png");
    std::cout << "Saved: plots/output/mass_ee_adaptive_sig_sim_ratios_coarse.{pdf,png}\n";

    // --------------------------------------------------------------------
    // Diagnostics: integrals (raw, undo /width)
    // --------------------------------------------------------------------
    auto rawIntegral = [](TH1D* h) {
        double s = 0.0;
        for (int b = 1; b <= h->GetNbinsX(); ++b)
            s += h->GetBinContent(b) * h->GetBinWidth(b);
        return s;
    };
    std::cout << "\n[sim] Weighted integrals (raw counts):\n"
              << "  none = " << rawIntegral(h_none) << "\n"
              << "  PT3  = " << rawIntegral(h_pt3)  << "\n"
              << "  PT2  = " << rawIntegral(h_pt2)  << "\n"
              << "  PT3/none = " << (rawIntegral(h_pt3) / rawIntegral(h_none)) << "\n"
              << "  PT2/PT3  = " << (rawIntegral(h_pt2) / rawIntegral(h_pt3))  << "\n";

    f->Close();
}
