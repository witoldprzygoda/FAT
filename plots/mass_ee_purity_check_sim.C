// mass_ee_purity_check_sim.C — research macro: how does the MC purity gate
// (ep_sim_id==2 && em_sim_id==3 && ep_sim_geninfo2==em_sim_geninfo2) reshape
// the reconstructed sim m_ee spectrum?
//
// WHY a separate research macro: the analysis output `output_epem_sim.root`
// has the purity gate HARD-CODED as an early return in main.cc:55-63, so
// nothing in the output represents the "without purity" reference. We have to
// read the raw SMASH source files directly via TChain (just like the earlier
// channel-discovery research macros).
//
// What is plotted (2×2 canvas, 1500 × 900 px):
//   Top-left  : REC m_ee spectrum, overlay of (no-purity) vs (with-purity)
//   Top-right : COR m_ee spectrum, overlay of (no-purity) vs (with-purity)
//   Bot-left  : ratio  N_REC(no-purity) / N_REC(with-purity)
//   Bot-right : ratio  N_COR(no-purity) / N_COR(with-purity)
//
// Quality cuts applied to BOTH versions (taken from main.cc, mirrors the
// real analysis selection apart from the purity gate):
//   isBest        == 1
//   eVertReco_z   >= -500 mm
//
// REC kinematics use (ep_p, ep_theta, ep_phi) and likewise for em — raw
// reconstructed magnitudes and angles, angles in DEGREES (HADES convention).
// COR kinematics use (ep_p_corr_ep, em_p_corr_em) — energy-loss corrected
// momenta under e+/e- hypothesis (same angles).  Both formulas use
// massive-electron 4-vector addition with m_e = 0.51099895 MeV/c² (see
// m_ee_GeV()).
//
// Yield-difference table (printed to stdout): per mass window (π⁰ /
// η Dalitz / high mass / full) plus REC/COR rows, showing weighted yields
// without vs with purity gate plus the fake fraction
//   fake_fraction = (Y_no − Y_with) / Y_no
// which quantifies how much of each m_ee region is "fake" (misidentified
// pairs or pairs from different decay vertices in multi-pair events).
//
// Default: chains 3 source files (≈132 M raw entries). Increase n_files
// up to 10 for the full sample (≈440 M entries → ~10 min runtime).
//
// Output:
//   plots/output/mass_ee_purity_check_sim.{pdf,png}
//
// Usage:
//   root -l -b -q plots/mass_ee_purity_check_sim.C
//   root -l -b -q 'plots/mass_ee_purity_check_sim.C(10)'   // full sample
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

// -------- m_ee histogram binning (matches the truth_vs_rec macro) --------
constexpr int    kNb   = 70;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.4;

// Physics constants for the m_ee calculation.
constexpr double kMeMeV = 0.51099895;
constexpr double kD2R   = 0.0174532925199432955;

// Quality cut thresholds, transcribed from setup_cuts.h (line 27-28):
//   cuts.defineValueCut("isBest", 1, "Best candidate selection");
//   cuts.defineMinCut("vertex_z", -500, "Vertex Z quality [mm]");
constexpr int    kIsBestVal = 1;
constexpr double kVertZMin  = -500.0;

// Colour code: "without purity" gets a warning red, "with purity" gets blue.
constexpr Color_t kColNoPur = kRed  + 1;
constexpr Color_t kColWPur  = kBlue + 1;

// Mass windows for the yield-difference table (same as the scan macros).
struct Window { std::string name; double lo; double hi; };
const std::vector<Window> kWindows = {
    {"A: pi0 (0-0.135)",        0.000, 0.135},
    {"B: eta Dal (0.135-0.6)",  0.135, 0.600},
    {"C: high mass (0.6-1.4)",  0.600, 1.400},
    {"D: full (0-1.4)",         0.000, 1.400},
};

// Compute m_ee from spherical (p [MeV], theta [deg], phi [deg]) for both legs
// using the full massive 4-vector formula. Returns m_ee in GeV/c².
double m_ee_GeV(double p1, double th1_deg, double ph1_deg,
                double p2, double th2_deg, double ph2_deg) {
    const double th1 = th1_deg * kD2R, ph1 = ph1_deg * kD2R;
    const double th2 = th2_deg * kD2R, ph2 = ph2_deg * kD2R;
    const double s1 = std::sin(th1), c1 = std::cos(th1);
    const double s2 = std::sin(th2), c2 = std::cos(th2);
    const double cphi = std::cos(ph1 - ph2);
    const double pdot = p1 * p2 * (s1 * s2 * cphi + c1 * c2);
    const double E1 = std::sqrt(p1 * p1 + kMeMeV * kMeMeV);
    const double E2 = std::sqrt(p2 * p2 + kMeMeV * kMeMeV);
    const double m2 = 2.0 * kMeMeV * kMeMeV + 2.0 * (E1 * E2 - pdot);
    if (m2 < 0.0) return 0.0;
    return std::sqrt(m2) * 1e-3;
}

std::vector<std::string> readList(const std::string& path, int n) {
    std::vector<std::string> out;
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "Cannot open list: " << path << "\n";
        return out;
    }
    std::string line;
    while (out.size() < static_cast<size_t>(n) && std::getline(f, line)) {
        line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
        line.erase(std::remove(line.begin(), line.end(), ','), line.end());
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back())))
            line.pop_back();
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.front())))
            line.erase(line.begin());
        if (!line.empty()) out.push_back(line);
    }
    return out;
}

void styleSpec(TH1D* h, Color_t color, Style_t marker) {
    h->SetLineColor(color);
    h->SetLineWidth(2);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.7);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.25);
}

// num / den with per-bin error propagation (same convention as elsewhere).
TH1D* makeRatio(TH1D* num, TH1D* den, const std::string& name) {
    TH1D* r = static_cast<TH1D*>(num->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    for (int b = 1; b <= num->GetNbinsX(); ++b) {
        const double n  = num->GetBinContent(b);
        const double en = num->GetBinError(b);
        const double d  = den->GetBinContent(b);
        const double ed = den->GetBinError(b);
        if (d != 0.0 && std::isfinite(d) && std::isfinite(n)) {
            const double val = n / d;
            const double rel = std::sqrt(
                (n != 0.0 ? std::pow(en / n, 2) : 0.0) +
                std::pow(ed / d, 2));
            r->SetBinContent(b, val);
            r->SetBinError  (b, std::abs(val) * rel);
        }
    }
    return r;
}

}  // anonymous namespace

void mass_ee_purity_check_sim(int n_files = 3) {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);

    // -------- Chain source files --------
    const std::string list_path =
        "/home/damian/hdd1/HADES/pp45/sim/smash/smash_lepton.list";
    auto files = readList(list_path, n_files);
    if (files.empty()) {
        std::cerr << "No SMASH source files found in " << list_path << "\n";
        return;
    }
    std::cout << "Chaining " << files.size()
              << " SMASH lepton source file(s):\n";
    TChain ch("EpEm_ID");
    for (const auto& p : files) {
        std::cout << "  " << p << "\n";
        ch.Add(p.c_str());
    }
    const Long64_t n_entries = ch.GetEntries();
    std::cout << "Total entries in chain: " << n_entries << "\n\n";

    // -------- TTreeReader values --------
    TTreeReader r(&ch);
    TTreeReaderValue<float> v_isBest(r, "isBest");
    TTreeReaderValue<float> v_vz    (r, "eVertReco_z");
    TTreeReaderValue<float> v_w     (r, "ep_sim_genweight");
    TTreeReaderValue<float> v_ep_id (r, "ep_sim_id");
    TTreeReaderValue<float> v_em_id (r, "em_sim_id");
    TTreeReaderValue<float> v_ep_g2 (r, "ep_sim_geninfo2");
    TTreeReaderValue<float> v_em_g2 (r, "em_sim_geninfo2");
    TTreeReaderValue<float> v_ep_p  (r, "ep_p");
    TTreeReaderValue<float> v_ep_th (r, "ep_theta");
    TTreeReaderValue<float> v_ep_ph (r, "ep_phi");
    TTreeReaderValue<float> v_em_p  (r, "em_p");
    TTreeReaderValue<float> v_em_th (r, "em_theta");
    TTreeReaderValue<float> v_em_ph (r, "em_phi");
    TTreeReaderValue<float> v_ep_p_cor(r, "ep_p_corr_ep");
    TTreeReaderValue<float> v_em_p_cor(r, "em_p_corr_em");

    // -------- Allocate histograms --------
    auto mkH = [](const char* n) {
        TH1D* h = new TH1D(n, "", kNb, kXmin, kXmax);
        h->Sumw2();
        h->SetDirectory(nullptr);
        return h;
    };
    TH1D* h_REC_no = mkH("h_REC_no");
    TH1D* h_REC_wp = mkH("h_REC_wp");
    TH1D* h_COR_no = mkH("h_COR_no");
    TH1D* h_COR_wp = mkH("h_COR_wp");

    // -------- Event loop --------
    Long64_t n_total = 0, n_qual = 0, n_pure = 0;
    while (r.Next()) {
        ++n_total;
        // Quality cuts (mirror main.cc's pre-purity gate)
        if (static_cast<int>(*v_isBest) != kIsBestVal) continue;
        if (*v_vz < kVertZMin)                          continue;
        ++n_qual;

        const double mee_REC = m_ee_GeV(*v_ep_p,    *v_ep_th, *v_ep_ph,
                                        *v_em_p,    *v_em_th, *v_em_ph);
        const double mee_COR = m_ee_GeV(*v_ep_p_cor,*v_ep_th, *v_ep_ph,
                                        *v_em_p_cor,*v_em_th, *v_em_ph);
        const double w = *v_w;

        // Always fill no-purity histograms.
        h_REC_no->Fill(mee_REC, w);
        h_COR_no->Fill(mee_COR, w);

        // Apply purity gate (sim_id + same-vertex).
        const bool pass_purity =
            (static_cast<int>(*v_ep_id) == 2)  &&
            (static_cast<int>(*v_em_id) == 3)  &&
            (*v_ep_g2 == *v_em_g2);
        if (!pass_purity) continue;
        ++n_pure;

        h_REC_wp->Fill(mee_REC, w);
        h_COR_wp->Fill(mee_COR, w);
    }

    std::cout << "Scanned:        " << n_total << "\n";
    std::cout << "Pass quality:   " << n_qual
              << "  ("
              << (100.0 * n_qual / std::max<Long64_t>(1, n_total))
              << " % of all)\n";
    std::cout << "Pass purity:    " << n_pure
              << "  ("
              << (100.0 * n_pure / std::max<Long64_t>(1, n_qual))
              << " % of qual)\n\n";

    // -------- Build ratios --------
    TH1D* r_REC = makeRatio(h_REC_no, h_REC_wp, "r_REC");
    TH1D* r_COR = makeRatio(h_COR_no, h_COR_wp, "r_COR");

    // -------- Per-window yield comparison --------
    std::cout << "Per-window yield comparison (weighted, sim_genweight):\n";
    printf("  %-26s   %-5s   %12s   %12s   %12s   %s\n",
           "window", "type", "no_purity", "with_purity",
           "diff (fake)", "fake_frac");
    std::cout << "  " << std::string(96, '-') << "\n";

    auto print_row = [](const char* win, const char* ver,
                        double y_no, double y_wp) {
        const double diff = y_no - y_wp;
        const double frac = (y_no > 0.0) ? diff / y_no * 100.0 : 0.0;
        printf("  %-26s   %-5s   %12.1f   %12.1f   %12.1f   %5.2f %%\n",
               win, ver, y_no, y_wp, diff, frac);
    };

    for (const auto& W : kWindows) {
        const int b_lo = h_REC_no->GetXaxis()->FindFixBin(W.lo + 1e-6);
        const int b_hi = h_REC_no->GetXaxis()->FindFixBin(W.hi - 1e-6);
        const double Y_REC_no = h_REC_no->Integral(b_lo, b_hi);
        const double Y_REC_wp = h_REC_wp->Integral(b_lo, b_hi);
        const double Y_COR_no = h_COR_no->Integral(b_lo, b_hi);
        const double Y_COR_wp = h_COR_wp->Integral(b_lo, b_hi);
        print_row(W.name.c_str(), "REC", Y_REC_no, Y_REC_wp);
        print_row("",             "COR", Y_COR_no, Y_COR_wp);
    }
    std::cout << "\n";

    // -------- Style --------
    styleSpec(h_REC_no, kColNoPur, 20);
    styleSpec(h_REC_wp, kColWPur,  21);
    styleSpec(h_COR_no, kColNoPur, 20);
    styleSpec(h_COR_wp, kColWPur,  21);
    styleSpec(r_REC,    kBlack, 20);
    styleSpec(r_COR,    kBlack, 20);

    // -------- Canvas: 2×2 --------
    TCanvas* c = new TCanvas("c_purity_check_sim",
        "Sim m_ee purity-gate impact (REC + COR)",
        1500, 900);
    c->Divide(2, 2, 0.001, 0.001);

    auto drawSpec = [&](int idx, TH1D* h_no, TH1D* h_wp,
                        const std::string& label) {
        c->cd(idx);
        gPad->SetLogy(true);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        const double y_max = std::max(h_no->GetMaximum(), h_wp->GetMaximum());
        h_no->SetTitle(Form(
            "m_{ee} (%s) — purity-gate impact;M_{e^{+}e^{-}} [GeV/c^{2}];weighted entries / 20 MeV",
            label.c_str()));
        h_no->GetYaxis()->SetRangeUser(std::max(1e-2, y_max * 1e-6),
                                       y_max * 5.0);
        h_no->Draw("E1");
        h_wp->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.45, 0.74, 0.95, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);
        leg->AddEntry(h_no,
            Form("no purity gate (#sum w = %.0f)", h_no->Integral()),  "lpe");
        leg->AddEntry(h_wp,
            Form("with sim_id + vertex (#sum w = %.0f)", h_wp->Integral()), "lpe");
        leg->Draw();
    };

    auto drawRatio = [&](int idx, TH1D* h, const std::string& label,
                         double y_lo, double y_hi) {
        c->cd(idx);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        h->SetTitle(Form(
            "Ratio %s_{no purity} / %s_{with purity};M_{e^{+}e^{-}} [GeV/c^{2}];N_{no} / N_{with}",
            label.c_str(), label.c_str()));
        h->GetYaxis()->SetRangeUser(y_lo, y_hi);
        h->Draw("E1");
        TLine* ln = new TLine(kXmin, 1.0, kXmax, 1.0);
        ln->SetLineStyle(2); ln->SetLineColor(kGray + 2); ln->SetLineWidth(1);
        ln->Draw();
    };

    drawSpec (1, h_REC_no, h_REC_wp, "REC");
    drawSpec (2, h_COR_no, h_COR_wp, "COR");
    // Range determined empirically: typical ratio sits at 1.000-1.005 in sim;
    // give some headroom in case there are noisier bins.
    drawRatio(3, r_REC, "REC", 0.95, 1.30);
    drawRatio(4, r_COR, "COR", 0.95, 1.30);

    gSystem->mkdir("plots/output", true);
    c->SaveAs("plots/output/mass_ee_purity_check_sim.pdf");
    c->SaveAs("plots/output/mass_ee_purity_check_sim.png");
    std::cout << "Saved: plots/output/mass_ee_purity_check_sim.{pdf,png}\n";
}
