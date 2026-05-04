// meson_dalitz_oa_vs_m.C — 2D scan of OA(e+e-) vs M(e+e-gamma).
//
// Reads from `meson_dalitz_nt` ntuple (mult==1, mass-window flags + RICH/ECAL
// diagnostics; filled WITHOUT the opening-angle cut so OA stays a free
// variable here).
//
//   X axis: M(e+e-gamma)  [GeV/c^2]   range [0.0, 1.2]   bin = 5 MeV
//   Y axis: OA(e+e-)      [deg]       range [0,   60]    bin = 0.5°
//
// Two outputs:
//   plots/output/meson_dalitz_oa_vs_m.{pdf,png}   — COLZ canvas (logz)
//   plots/output/meson_dalitz_oa_vs_m.root        — TH2D for projections/slices
//
// Slice example (interactive ROOT):
//   TFile f("plots/output/meson_dalitz_oa_vs_m.root");
//   auto h = (TH2D*)f.Get("h_oa_vs_m_epemg");
//   auto px = h->ProjectionX("px_pi0", h->GetYaxis()->FindBin(4.0),
//                                       h->GetYaxis()->FindBin(10.0));
//
// Usage: root -l -b -q plots/meson_dalitz_oa_vs_m.C

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>

void meson_dalitz_oa_vs_m() {

    TFile* f = TFile::Open("output_epem_exp.root", "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open output_epem_exp.root\n";
        return;
    }
    auto* t = (TTree*)f->Get("meson_dalitz_nt");
    if (!t) {
        std::cerr << "No meson_dalitz_nt ntuple in output_epem_exp.root\n";
        f->Close();
        return;
    }

    // Binning
    const int    nx = 240;     // 5 MeV/bin
    const double xmin = 0.0, xmax = 1.2;
    const int    ny = 120;     // 0.5°/bin
    const double ymin = 0.0, ymax = 60.0;

    TH2D* h = new TH2D("h_oa_vs_m_epemg",
        "#theta_{open}(e^{+}e^{-}) vs M(e^{+}e^{-}#gamma)  "
        "(mult==1, after basic event cuts);"
        "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];"
        "#theta_{open}(e^{+}e^{-}) [deg]",
        nx, xmin, xmax, ny, ymin, ymax);
    h->Sumw2();
    // NOTE: do NOT call SetDirectory(nullptr) here — TTree::Draw with the
    // "Y:X>>name" syntax looks the histogram up in gDirectory by name. We
    // detach AFTER the fill, just before the input file is closed.

    // Fill from ntuple. TTree::Draw uses "Y:X" convention.
    t->Draw("oa_epem:m_epemg>>h_oa_vs_m_epemg", "", "goff");
    h->SetDirectory(nullptr);  // detach so it survives f->Close() below
    std::cout << "Filled 2D histogram: " << h->GetEntries() << " entries\n";

    // Style + COLZ canvas with log Z.
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(100);

    auto* c = new TCanvas("c_meson_dalitz_oa_vs_m", "OA vs M(epemg)", 900, 700);
    c->SetMargin(0.12, 0.14, 0.12, 0.08);
    c->SetLogz();

    h->Draw("COLZ");
    c->Update();

    gSystem->mkdir("plots/output", kTRUE);
    c->SaveAs("plots/output/meson_dalitz_oa_vs_m.pdf");
    c->SaveAs("plots/output/meson_dalitz_oa_vs_m.png");

    // Persist the 2D histogram so it can be reopened for projections
    // and 1D slice analysis without re-running the macro on the ntuple.
    TFile* fout = TFile::Open("plots/output/meson_dalitz_oa_vs_m.root", "RECREATE");
    h->Write();
    fout->Close();

    std::cout << "Saved: plots/output/meson_dalitz_oa_vs_m.{pdf,png,root}\n";

    f->Close();
}
