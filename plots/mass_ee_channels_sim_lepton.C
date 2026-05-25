// mass_ee_channels_sim_lepton.C — RESEARCH macro for SMASH LEPTON channel
// decomposition. Reads SMASH source files DIRECTLY (output_epem_sim.root does
// not carry sim_parentid / sim_geninfo* branches), classifies each lepton pair
// by its production channel via `sim_geninfo1`, and plots m_ee_sim broken
// down by channel.
//
// Source: /home/damian/hdd1/HADES/pp45/sim/smash/smash_lepton.list
// Tree:   EpEm_ID
//
// Channel encoding observed in EpEm_ID (sim_geninfo1):
//   7051  (~78 %) — π⁰ Dalitz   (π⁰ → e⁺e⁻γ)
//   17051 (~9.5%) — η  Dalitz   (η  → e⁺e⁻γ)
//   41    (~7.5%) — ω           (ω  → e⁺e⁻ direct OR ω → π⁰ e⁺e⁻ Dalitz)
//   52    (~3.7%) — likely ρ / Δ (TBD — examine m_ee shape)
//   7     (~0.8%) — π⁰ γγ → e⁺e⁻ conversion (γ parent, π⁰ grandparent)
//   55    (~0.2%) — TBD
//   17    (~0.1%) — η  γγ → e⁺e⁻ conversion (γ parent, η  grandparent)
//
// Each lepton is essentially primary in SMASH (parentid/grandparentid usually
// 0 because the meson decay happens in PLUTO, not in GEANT cascade). The
// physics origin is encoded in `sim_geninfo1`. processid==1 = generator
// primary, ==6 = pair-production conversion in detector material, ==7 = other
// secondary process.
//
// Three panels:
//   Pad 1 — geninfo1 distribution (sorted, weighted) on log Y
//   Pad 2 — m_ee_sim overlay (one curve per major channel + total) on log Y
//   Pad 3 — m_ee_sim stack (THStack) on log Y
//
// Arguments:
//   1) n_files     — number of source files to chain (default 1; max 10)
//   2) use_weights — true → fillw(sim_genweight); false → fill(1.0); default true
//
// Usage:
//   root -l -b -q plots/mass_ee_channels_sim_lepton.C
//   root -l -b -q 'plots/mass_ee_channels_sim_lepton.C(3)'
//   root -l -b -q 'plots/mass_ee_channels_sim_lepton.C(10, false)'
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

constexpr double kMeMeV = 0.51099895;    // electron mass in MeV/c²
constexpr int    kNb    = 200;
constexpr double kXmin  = 0.0;
constexpr double kXmax  = 1.4;           // m_ee_sim in GeV/c²

// Channel definition: name, geninfo1 value (or -1 = catch-all), color, marker.
struct Channel {
    std::string name;
    int         geninfo1;        // exact match (or -1 for "other")
    Color_t     color;
    Style_t     marker;
};

// Order = legend order (top-down). The "other" entry catches anything that
// did not match an earlier rule (its geninfo1 is ignored).
// Encoding pattern observed in the data:
//   X * 1000 + 51   → meson X Dalitz decay (3-body: X → e⁺e⁻γ or X → e⁺e⁻N)
//   X (small)       → meson X direct decay to e⁺e⁻ (2-body)
//   7 / 17 (small)  → π⁰ / η  γγ → e⁺e⁻ conversion in detector material
// Guesses for X (HADES/HGeant convention): 7=π⁰, 17=η, 41=ω, 52=?(ρ⁰ or Δ),
// 55=?. Confirm or correct with the m_ee shape.
const std::vector<Channel> kChannels = {
    {"#pi^{0} Dalitz (7051)",         7051, kRed     + 1, 20},
    {"#pi^{0} #gamma#gamma conv. (7)",   7, kRed     + 3, 24},
    {"#eta Dalitz (17051)",          17051, kBlue    + 1, 21},
    {"#eta #gamma#gamma conv. (17)",    17, kBlue    + 3, 25},
    {"ch 52 Dalitz (52051)",         52051, kMagenta + 2, 23},
    {"ch 41 direct e^{+}e^{-}",         41, kGreen   + 2, 22},
    {"ch 52 direct e^{+}e^{-}",         52, kMagenta + 1, 32},
    {"ch 55 direct e^{+}e^{-}",         55, kOrange  + 7, 33},
    {"other",                           -1, kGray    + 2, 34},
};

double m_ee_sim_GeV(double pxa, double pya, double pza,
                    double pxb, double pyb, double pzb) {
    const double E1 = std::sqrt(pxa*pxa + pya*pya + pza*pza + kMeMeV*kMeMeV);
    const double E2 = std::sqrt(pxb*pxb + pyb*pyb + pzb*pzb + kMeMeV*kMeMeV);
    const double dot = pxa*pxb + pya*pyb + pza*pzb;
    const double m2  = 2.0*kMeMeV*kMeMeV + 2.0*(E1*E2 - dot);
    if (m2 < 0.0) return 0.0;
    return std::sqrt(m2) * 1e-3;            // MeV → GeV
}

int classify(int geninfo1) {
    for (size_t i = 0; i < kChannels.size(); ++i) {
        if (kChannels[i].geninfo1 == geninfo1) return static_cast<int>(i);
    }
    // "other" bucket is the last entry whose geninfo1 == -1
    for (size_t i = 0; i < kChannels.size(); ++i)
        if (kChannels[i].geninfo1 == -1) return static_cast<int>(i);
    return -1;
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

void styleHist(TH1D* h, Color_t c, Style_t m) {
    h->SetLineColor(c);
    h->SetLineWidth(2);
    h->SetMarkerColor(c);
    h->SetMarkerStyle(m);
    h->SetMarkerSize(0.7);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.040);
    h->GetYaxis()->SetLabelSize(0.040);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.15);
}

}  // anonymous namespace

void mass_ee_channels_sim_lepton(int n_files = 1, bool use_weights = true) {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    const std::string list_path =
        "/home/damian/hdd1/HADES/pp45/sim/smash/smash_lepton.list";
    auto files = readList(list_path, n_files);
    if (files.empty()) {
        std::cerr << "No input files found\n";
        return;
    }
    std::cout << "Chaining " << files.size() << " SMASH lepton source file(s):\n";
    TChain ch("EpEm_ID");
    for (const auto& p : files) {
        std::cout << "  " << p << "\n";
        ch.Add(p.c_str());
    }
    const Long64_t n_entries = ch.GetEntries();
    std::cout << "Total entries in chain: " << n_entries << "\n";

    TTreeReader r(&ch);
    // PID purity gate
    TTreeReaderValue<float> ep_sim_id  (r, "ep_sim_id");
    TTreeReaderValue<float> em_sim_id  (r, "em_sim_id");
    TTreeReaderValue<float> ep_g2      (r, "ep_sim_geninfo2");
    TTreeReaderValue<float> em_g2      (r, "em_sim_geninfo2");
    // Channel classifier — use ep side (matched at gate by g2 equality)
    TTreeReaderValue<float> ep_g1      (r, "ep_sim_geninfo1");
    // Weight
    TTreeReaderValue<float> ep_w       (r, "ep_sim_genweight");
    // SIM momenta (MeV)
    TTreeReaderValue<float> ep_px      (r, "ep_sim_px");
    TTreeReaderValue<float> ep_py      (r, "ep_sim_py");
    TTreeReaderValue<float> ep_pz      (r, "ep_sim_pz");
    TTreeReaderValue<float> em_px      (r, "em_sim_px");
    TTreeReaderValue<float> em_py      (r, "em_sim_py");
    TTreeReaderValue<float> em_pz      (r, "em_sim_pz");

    // Allocate histograms
    const int nch = static_cast<int>(kChannels.size());
    std::vector<TH1D*> h_mee(nch, nullptr);
    for (int i = 0; i < nch; ++i) {
        h_mee[i] = new TH1D(Form("h_mee_ch%d", i),
                            "m_{e^{+}e^{-}}^{sim} per channel;"
                            "M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];"
                            "weighted entries / 7 MeV",
                            kNb, kXmin, kXmax);
        h_mee[i]->Sumw2();
        styleHist(h_mee[i], kChannels[i].color, kChannels[i].marker);
    }
    TH1D* h_mee_tot = new TH1D("h_mee_tot",
                               "m_{e^{+}e^{-}}^{sim} total;"
                               "M_{e^{+}e^{-}}^{sim} [GeV/c^{2}];"
                               "weighted entries / 7 MeV",
                               kNb, kXmin, kXmax);
    h_mee_tot->Sumw2();
    styleHist(h_mee_tot, kBlack, 20);
    h_mee_tot->SetLineWidth(3);

    // Channel counters
    std::vector<double> n_w(nch, 0.0), n_u(nch, 0.0);
    double n_w_tot = 0.0, n_u_tot = 0.0;

    // Also keep an unbinned tally of geninfo1 values for the top-N table.
    std::map<int, std::pair<double,double>> by_g1; // geninfo1 -> (sum_w, count)

    Long64_t pass = 0, scanned = 0;
    while (r.Next()) {
        ++scanned;
        if (static_cast<int>(*ep_sim_id) != 2)            continue;
        if (static_cast<int>(*em_sim_id) != 3)            continue;
        if (*ep_g2 != *em_g2)                              continue;
        ++pass;

        const double w = use_weights ? static_cast<double>(*ep_w) : 1.0;
        const int    g1 = static_cast<int>(*ep_g1);
        const int    ci = classify(g1);
        const double m  = m_ee_sim_GeV(*ep_px, *ep_py, *ep_pz,
                                       *em_px, *em_py, *em_pz);

        h_mee_tot->Fill(m, w);
        if (ci >= 0) {
            h_mee[ci]->Fill(m, w);
            n_w[ci] += w;
            n_u[ci] += 1.0;
        }
        n_w_tot += w;
        n_u_tot += 1.0;

        auto& slot = by_g1[g1];
        slot.first  += w;
        slot.second += 1.0;
    }
    std::cout << "Scanned: " << scanned << "   Passed purity+vertex: "
              << pass << "  (" << (100.0 * pass / std::max<Long64_t>(1, scanned))
              << " %)\n";
    std::cout << "Weight mode: " << (use_weights ? "WEIGHTED (sim_genweight)" : "UNWEIGHTED (1.0)") << "\n\n";

    // --- Table: channel breakdown (named) ---
    std::cout << "Per-channel breakdown (weighted, then unweighted counts):\n";
    std::cout << "  " << std::string(64, '-') << "\n";
    std::cout << "  " << "channel                              "
              << "    sum_w       % w    " << "   count       % N\n";
    std::cout << "  " << std::string(64, '-') << "\n";
    for (int i = 0; i < nch; ++i) {
        printf("  %-38s  %10.1f  %6.2f%%   %10.0f  %6.2f%%\n",
               kChannels[i].name.c_str(),
               n_w[i], 100.0 * n_w[i] / std::max(1e-9, n_w_tot),
               n_u[i], 100.0 * n_u[i] / std::max(1e-9, n_u_tot));
    }
    std::cout << "  " << std::string(64, '-') << "\n";
    printf("  %-38s  %10.1f  %6.2f%%   %10.0f  %6.2f%%\n",
           "TOTAL", n_w_tot, 100.0, n_u_tot, 100.0);
    std::cout << "\n";

    // --- Table: top distinct geninfo1 values (raw) ---
    std::vector<std::pair<int, std::pair<double,double>>> sorted(
        by_g1.begin(), by_g1.end());
    std::sort(sorted.begin(), sorted.end(),
              [](const auto& a, const auto& b) {
                  return a.second.first > b.second.first;
              });
    std::cout << "Top 12 distinct sim_geninfo1 values (sorted by weighted sum):\n";
    std::cout << "  " << std::string(48, '-') << "\n";
    std::cout << "  geninfo1      sum_w       %       count\n";
    std::cout << "  " << std::string(48, '-') << "\n";
    for (size_t i = 0; i < sorted.size() && i < 12; ++i) {
        printf("  %8d  %12.1f  %6.2f%%  %10.0f\n",
               sorted[i].first,
               sorted[i].second.first,
               100.0 * sorted[i].second.first / std::max(1e-9, n_w_tot),
               sorted[i].second.second);
    }
    std::cout << "\n";

    // ====================================================================
    // Plot 1 — Channel bar chart (top N geninfo1)
    // ====================================================================
    const int top_n = std::min<int>(10, sorted.size());
    TH1D* h_bar = new TH1D("h_bar",
                           "Top sim_geninfo1 channels (weighted);"
                           "sim_geninfo1;sum of sim_genweight",
                           top_n, 0, top_n);
    for (int i = 0; i < top_n; ++i) {
        h_bar->SetBinContent(i + 1, sorted[i].second.first);
        h_bar->GetXaxis()->SetBinLabel(
            i + 1, Form("%d", sorted[i].first));
    }
    h_bar->SetFillColor(kAzure - 4);
    h_bar->SetLineColor(kAzure - 7);
    styleHist(h_bar, kAzure - 7, 1);
    h_bar->SetMarkerStyle(0);

    // ====================================================================
    // Canvas: 2 pads — bar chart + m_ee overlay (log Y)
    // ====================================================================
    TCanvas* c = new TCanvas(
        "c_channels",
        Form("Channel decomposition (n_files=%d, %s)",
             (int)files.size(), use_weights ? "weighted" : "unweighted"),
        1600, 700);
    c->Divide(2, 1, 0.001, 0.001);

    // ---- Pad 1: bar chart of top geninfo1 ----
    c->cd(1);
    gPad->SetMargin(0.13, 0.04, 0.16, 0.10);
    gPad->SetLogy(true);
    h_bar->GetYaxis()->SetRangeUser(
        std::max(1.0, sorted.back().second.first * 0.5),
        sorted.front().second.first * 4.0);
    h_bar->GetXaxis()->LabelsOption("h");
    h_bar->Draw("BAR HIST");
    TLatex lt;
    lt.SetNDC(); lt.SetTextSize(0.040);
    lt.DrawLatex(0.16, 0.92,
                 Form("SMASH lepton — %lld pairs after purity+vertex gate",
                      static_cast<long long>(pass)));

    // ---- Pad 2: m_ee_sim overlay by channel ----
    c->cd(2);
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
    gPad->SetLogy(true);

    // Y range: from total max to a sensible floor.
    const double y_max = h_mee_tot->GetMaximum();
    const double y_min = std::max(1e-3, y_max * 1e-6);
    h_mee_tot->GetYaxis()->SetRangeUser(y_min, y_max * 5.0);
    h_mee_tot->Draw("HIST");
    for (int i = 0; i < nch; ++i) {
        if (h_mee[i]->GetEntries() == 0) continue;
        h_mee[i]->Draw("HIST SAME");
    }

    TLegend* leg = new TLegend(0.55, 0.55, 0.95, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.028);
    leg->AddEntry(h_mee_tot, "total", "l");
    for (int i = 0; i < nch; ++i) {
        if (h_mee[i]->GetEntries() == 0) continue;
        leg->AddEntry(h_mee[i],
                      Form("%s  (%.1f%%)",
                           kChannels[i].name.c_str(),
                           100.0 * n_w[i] / std::max(1e-9, n_w_tot)),
                      "l");
    }
    leg->Draw();

    // -------- Save --------
    gSystem->mkdir("plots/output", true);
    const std::string nf_suffix = Form("_n%d", (int)files.size());
    const std::string wt_suffix = use_weights ? std::string{} : "_unw";
    const std::string base =
        std::string("plots/output/mass_ee_channels_sim_lepton") +
        nf_suffix + wt_suffix;
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
    std::cout << "Saved: " << base << ".{pdf,png}\n";
}
