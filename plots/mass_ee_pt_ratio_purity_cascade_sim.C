// mass_ee_pt_ratio_purity_cascade_sim.C
// =========================================================================
// Cumulative cut cascade study on sim m_ee spectrum and PT2/PT3 trigger ratio.
// Five scenarios, each successively tighter:
//
//   1) RAW                — all reconstructed e⁺e⁻ candidates, no cuts
//   2) +eVertReco_z > -500
//   3) +isBest == 1
//   4) +ep_sim_id == 2 && em_sim_id == 3   (Geant3 PID purity)
//   5) +epem_same_vertex == 1               (= ep_sim_geninfo2 == em_sim_geninfo2)
//
// Step 5 reproduces the previously hard-gated sample (286M events).  Steps 1-4
// expose progressively contaminated samples (fake pairs, misidentified hadrons,
// random-vertex combinations) that were previously thrown away at fill time.
//
// All histograms are weighted by sim_genweight (per-event SMASH luminosity).
//
// Per scenario, one 1×2 canvas:
//   Left   : m_ee spectrum, PT3 (black filled circle) overlaid with PT2 (red
//            filled square), log Y.
//   Right  : ratio N_PT2 / N_PT3 with per-bin error propagation, constant fit
//            in [0.1, 1.1] GeV/c² drawn in blue. Y axis fixed [0, 3] for
//            cross-scenario comparability. Reference y=1 dashed grey line.
//
// At the end: stdout summary table of (a, chi²/ndf, PT3 yield, PT2 yield).
//
// Output:
//   plots/output/mass_ee_pt_ratio_purity_cascade_sim_step1_raw.{pdf,png}
//   plots/output/mass_ee_pt_ratio_purity_cascade_sim_step2_vertex.{pdf,png}
//   plots/output/mass_ee_pt_ratio_purity_cascade_sim_step3_isBest.{pdf,png}
//   plots/output/mass_ee_pt_ratio_purity_cascade_sim_step4_simID.{pdf,png}
//   plots/output/mass_ee_pt_ratio_purity_cascade_sim_step5_samevtx.{pdf,png}
//   plots/output/mass_ee_pt_ratio_purity_cascade_sim_overlay.{pdf,png}
//     (extra: overlay of all 5 ratios on one panel for direct comparison)
//
// Usage:
//   root -l -b -q plots/mass_ee_pt_ratio_purity_cascade_sim.C
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

constexpr const char* kInputFile = "output_epem_sim.root";
constexpr const char* kNtName    = "dilepton_nt";         // m_ee = m_ee_rec here

constexpr int    kNb    = 70;
constexpr double kXmin  = 0.0;
constexpr double kXmax  = 1.4;
constexpr double kFitLo = 0.1;
constexpr double kFitHi = 1.1;

constexpr Color_t kColPT3   = kBlack;
constexpr Color_t kColPT2   = kRed + 1;
constexpr Color_t kColRatio = kBlue + 1;
constexpr Color_t kColFit   = kBlue + 2;

struct Scenario {
    std::string id;
    std::string label;
    std::string cut;       // applied as MULTIPLICATIVE factor in the weight string
    std::string tree;      // tree name to read from (default: dilepton_nt)
    std::string mass_var;  // mass variable to draw (default: m_ee)
};

const std::vector<Scenario> kScenarios = {
    {"step0_raw_truth", "0) RAW truth — no cuts",                          "1", "dilepton_nt_cor", "m_ee_sim"},
    {"step1_raw",       "1) RAW — no cuts",                                "1", "dilepton_nt", "m_ee"},
    {"step2_vertex",    "2) +eVertReco_z > -500",                          "(eVertReco_z>-500)", "dilepton_nt", "m_ee"},
    {"step3_isBest",    "3) +isBest == 1",                                  "(eVertReco_z>-500)*(isBest==1)", "dilepton_nt", "m_ee"},
    {"step4_simID",     "4) +ep_sim_id==2 && em_sim_id==3",                "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)", "dilepton_nt", "m_ee"},
    {"step5_samevtx",   "5) +epem_same_vertex==1 (full old purity gate)",  "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)*(epem_same_vertex==1)", "dilepton_nt", "m_ee"},
};

// Convention colors for the OVERLAY (one per scenario; gradient for clarity).
// step0 (truth, raw) is kGray+2 so it is visually distinct from the rec-based scenarios.
const std::array<Color_t, 6> kOverlayColors = {
    kGray + 2, kRed + 1, kOrange + 7, kGreen + 2, kAzure + 1, kViolet + 2
};

TH1D* drawWithWeight(TTree* t, const std::string& mass_var,
                     const std::string& weight, const std::string& name) {
    // Create histogram IN gDirectory so TTree::Draw fills it (not a new one).
    TH1D* h = new TH1D(name.c_str(), "", kNb, kXmin, kXmax);
    h->Sumw2();
    t->Draw((mass_var + ">>" + name).c_str(),
            weight.c_str(), "goff");
    h->SetDirectory(nullptr);     // detach AFTER fill so it survives file close
    return h;
}

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

void styleHist(TH1D* h, Color_t c, Style_t m) {
    h->SetLineColor(c);
    h->SetMarkerColor(c);
    h->SetMarkerStyle(m);
    h->SetMarkerSize(0.8);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

struct FitResult { double a; double e; double chi2_ndf; int ndf; };

FitResult fitConst(TH1D* h, double xlo, double xhi, Color_t color) {
    const std::string fname = std::string(h->GetName()) + "_fc";
    TF1* f = new TF1(fname.c_str(), "[0]", xlo, xhi);
    f->SetLineColor(color);
    f->SetLineWidth(3);
    h->Fit(f, "RQ");          // R = use range, Q = quiet (attaches the fit)
    FitResult r;
    r.a        = f->GetParameter(0);
    r.e        = f->GetParError(0);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    return r;
}

}  // anonymous namespace

void mass_ee_pt_ratio_purity_cascade_sim() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);

    TFile* fs = TFile::Open(kInputFile, "READ");
    if (!fs || fs->IsZombie()) {
        std::cerr << "Cannot open " << kInputFile << "\n"; return;
    }
    TTree* nt = dynamic_cast<TTree*>(fs->Get(kNtName));
    if (!nt) {
        std::cerr << kNtName << " missing in " << kInputFile << "\n"; return;
    }
    // Cache of tree pointers keyed by name. dilepton_nt_cor is loaded on demand
    // (only scenario 0 needs it; missing tree would be fatal but useful to know).
    std::map<std::string, TTree*> tree_cache;
    tree_cache[kNtName] = nt;
    {
        TTree* nt_cor = dynamic_cast<TTree*>(fs->Get("dilepton_nt_cor"));
        if (!nt_cor) {
            std::cerr << "dilepton_nt_cor missing in " << kInputFile << "\n"; return;
        }
        tree_cache["dilepton_nt_cor"] = nt_cor;
    }
    std::cout << "Input: " << kInputFile << " / " << kNtName
              << " (" << nt->GetEntries() << " entries)\n";
    std::cout << "Aux:   " << kInputFile << " / dilepton_nt_cor"
              << " (" << tree_cache["dilepton_nt_cor"]->GetEntries() << " entries)\n";
    std::cout << "Fit range: [" << kFitLo << ", " << kFitHi << "] GeV/c²\n";
    std::cout << "Weight: sim_genweight × cumulative cut\n\n";

    gSystem->mkdir("plots/output", true);

    // Storage for the overlay canvas + summary table.
    std::vector<TH1D*> ratios_for_overlay;
    std::vector<FitResult> fits;
    std::vector<double> yields_PT3, yields_PT2;

    for (size_t i = 0; i < kScenarios.size(); ++i) {
        const auto& S = kScenarios[i];
        std::cout << "\n=== Scenario " << S.label << " ===\n";

        const std::string w_PT3 = "(pt3==1)*" + S.cut + "*sim_genweight";
        const std::string w_PT2 = "(pt2==1)*" + S.cut + "*sim_genweight";

        std::cout << "  tree: " << S.tree << "   mass var: " << S.mass_var << "\n";
        std::cout << "  PT3 weight: " << w_PT3 << "\n";

        // Pick the appropriate tree for this scenario from the cache.
        TTree* t = tree_cache[S.tree];

        // -------- Fill PT3 + PT2 histograms with cumulative weight --------
        TH1D* h_PT3 = drawWithWeight(t, S.mass_var, w_PT3,
            Form("h_PT3_%s", S.id.c_str()));
        TH1D* h_PT2 = drawWithWeight(t, S.mass_var, w_PT2,
            Form("h_PT2_%s", S.id.c_str()));

        const double Y3 = h_PT3->Integral();
        const double Y2 = h_PT2->Integral();
        std::cout << "  yields: PT3 = " << Y3 << "  PT2 = " << Y2 << "\n";
        yields_PT3.push_back(Y3);
        yields_PT2.push_back(Y2);

        // -------- Ratio + constant fit --------
        TH1D* r = makeRatio(h_PT2, h_PT3,
            Form("r_%s", S.id.c_str()));

        styleHist(h_PT3, kColPT3, 20);     // filled circle
        styleHist(h_PT2, kColPT2, 21);     // filled square
        styleHist(r,     kColRatio, 20);

        FitResult fr = fitConst(r, kFitLo, kFitHi, kColFit);
        printf("  fit [%.2f, %.2f]:  a = %.4f ± %.4f   chi2/ndf = %.1f / %d = %.2f\n",
               kFitLo, kFitHi, fr.a, fr.e, fr.chi2_ndf * fr.ndf, fr.ndf, fr.chi2_ndf);
        fits.push_back(fr);

        // -------- Canvas: 1×2 --------
        TCanvas* c = new TCanvas(
            Form("c_purity_cascade_%s", S.id.c_str()),
            S.label.c_str(),
            1500, 600);
        c->Divide(2, 1, 0.001, 0.001);

        // Left: m_ee log Y
        c->cd(1);
        gPad->SetLogy(true);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        const double y_max = std::max(h_PT3->GetMaximum(), h_PT2->GetMaximum());
        h_PT3->SetTitle(Form(
            "m_{ee} spectrum — %s;M_{e^{+}e^{-}} [GeV/c^{2}];weighted entries / 20 MeV",
            S.label.c_str()));
        h_PT3->GetYaxis()->SetRangeUser(std::max(1e-2, y_max * 1e-6), y_max * 5.0);
        h_PT3->Draw("E1");
        h_PT2->Draw("E1 SAME");
        TLegend* leg = new TLegend(0.55, 0.74, 0.95, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
        leg->AddEntry(h_PT3, Form("PT3  (#sum w = %.0f)", Y3), "lpe");
        leg->AddEntry(h_PT2, Form("PT2  (#sum w = %.0f)", Y2), "lpe");
        leg->Draw();

        // Right: ratio with constant fit
        c->cd(2);
        gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
        r->SetTitle(Form(
            "N_{PT2}/N_{PT3} — %s;M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT2}/N_{PT3}",
            S.label.c_str()));
        r->GetYaxis()->SetRangeUser(0.5, 2.0);
        r->Draw("E1");          // fit drawn automatically (TF1 attached)
        TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
        lref->SetLineStyle(2); lref->SetLineColor(kGray + 2); lref->Draw();

        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.034);
        tx.SetTextColor(kColFit);
        tx.DrawLatex(0.16, 0.86,
            Form("const fit [%.2f, %.2f]:  a = %.4f #pm %.4f",
                 kFitLo, kFitHi, fr.a, fr.e));
        tx.SetTextSize(0.030);
        tx.SetTextColor(kGray + 3);
        tx.DrawLatex(0.16, 0.81,
            Form("#chi^{2}/ndf = %.1f / %d = %.2f",
                 fr.chi2_ndf * fr.ndf, fr.ndf, fr.chi2_ndf));

        // ----- Save -----
        const std::string base =
            std::string("plots/output/mass_ee_pt_ratio_purity_cascade_sim_") + S.id;
        c->SaveAs((base + ".pdf").c_str());
        c->SaveAs((base + ".png").c_str());
        std::cout << "  saved: " << base << ".{pdf,png}\n";

        // Store ratio for overlay canvas (clone so styling doesn't conflict).
        TH1D* r_ovl = static_cast<TH1D*>(r->Clone(Form("r_ovl_%s", S.id.c_str())));
        r_ovl->SetDirectory(nullptr);
        // Strip the attached fit so the overlay marker series is clean.
        r_ovl->GetListOfFunctions()->Clear();
        styleHist(r_ovl, kOverlayColors[i], 20 + (int)i);
        ratios_for_overlay.push_back(r_ovl);
    }

    // ============================================================================
    // Overlay canvas — all 5 ratios on one panel + summary table
    // ============================================================================
    TCanvas* c_ovl = new TCanvas("c_purity_cascade_overlay",
        "All 5 ratios overlaid", 1500, 700);
    gPad->SetMargin(0.10, 0.04, 0.13, 0.10);
    ratios_for_overlay[0]->SetTitle(
        "PT2/PT3 cumulative purity cascade — sim;"
        "M_{e^{+}e^{-}} [GeV/c^{2}];N_{PT2}/N_{PT3}");
    ratios_for_overlay[0]->GetYaxis()->SetRangeUser(0.5, 2.5);
    ratios_for_overlay[0]->Draw("E1");
    for (size_t i = 1; i < ratios_for_overlay.size(); ++i)
        ratios_for_overlay[i]->Draw("E1 SAME");

    TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
    lref->SetLineStyle(2); lref->SetLineColor(kGray + 2); lref->Draw();

    TLegend* leg = new TLegend(0.55, 0.62, 0.96, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.022);
    for (size_t i = 0; i < kScenarios.size(); ++i) {
        leg->AddEntry(ratios_for_overlay[i],
            Form("%s   a=%.3f #chi^{2}/ndf=%.1f",
                 kScenarios[i].label.c_str(),
                 fits[i].a, fits[i].chi2_ndf),
            "lpe");
    }
    leg->Draw();

    c_ovl->SaveAs("plots/output/mass_ee_pt_ratio_purity_cascade_sim_overlay.pdf");
    c_ovl->SaveAs("plots/output/mass_ee_pt_ratio_purity_cascade_sim_overlay.png");

    // ============================================================================
    // Summary table
    // ============================================================================
    std::cout << "\n\n=== SUMMARY ===\n";
    printf("  %-50s  %-12s  %-12s  %-12s  %-15s\n",
           "scenario", "PT3 yield", "PT2 yield", "PT2/PT3 a", "chi2/ndf");
    std::cout << "  " << std::string(110, '-') << "\n";
    for (size_t i = 0; i < kScenarios.size(); ++i) {
        printf("  %-50s  %12.0f  %12.0f  %.4f±%.4f  %.1f / %d = %.2f\n",
               kScenarios[i].label.c_str(),
               yields_PT3[i], yields_PT2[i],
               fits[i].a, fits[i].e,
               fits[i].chi2_ndf * fits[i].ndf, fits[i].ndf, fits[i].chi2_ndf);
    }
    std::cout << "\n";

    fs->Close();
}
