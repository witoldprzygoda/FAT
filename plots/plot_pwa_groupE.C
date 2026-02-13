/**
 * @file plot_pwa_groupE.C
 * @brief Plot PWA Group E: Gottfried-Jackson angle distributions
 *
 * Displays pwa_pip_gj, pwa_pipn_gj, pwa_n_gj from three files
 * File 1 (black): output_fwd_ppip.root
 * File 2 (blue): output_fwd_ppip_full.root
 * File 3 (red): output_fwd_ppip_neutroncut.root
 */

void plot_pwa_groupE() {
    gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.12);

    // Open input files
    TFile* f1 = TFile::Open("output_fwd_ppip.root");
    TFile* f2 = TFile::Open("output_fwd_ppip_full.root");
    TFile* f3 = TFile::Open("output_fwd_ppip_neutroncut.root");

    if (!f1 || !f2 || !f3) {
        std::cerr << "Error: Could not open input files!" << std::endl;
        return;
    }

    // Histogram names for Group E
    const char* hist_names[3] = {"pwa/pwa_pip_gj", "pwa/pwa_pipn_gj", "pwa/pwa_n_gj"};
    const char* titles[3] = {"#pi^{+} GJ angle", "#pi^{+}n GJ angle", "Neutron GJ angle"};

    // Create canvas with 3 horizontal pads
    TCanvas* c1 = new TCanvas("c1", "PWA Group E: Gottfried-Jackson angle distributions", 1800, 600);
    c1->Divide(3, 1);

    // Loop over three histograms
    for (int i = 0; i < 3; ++i) {
        c1->cd(i + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);

        // Get histograms from three files
        TH1D* h1 = (TH1D*)f1->Get(hist_names[i]);
        TH1D* h2 = (TH1D*)f2->Get(hist_names[i]);
        TH1D* h3 = (TH1D*)f3->Get(hist_names[i]);

        if (!h1 || !h2 || !h3) {
            std::cerr << "Error: Could not get histogram " << hist_names[i] << std::endl;
            continue;
        }

        // Set colors and markers
        h1->SetLineColor(kBlack);
        h1->SetMarkerColor(kBlack);
        h1->SetMarkerStyle(20);
        h1->SetMarkerSize(0.8);
        h1->SetLineWidth(1);

        h2->SetLineColor(kBlue);
        h2->SetMarkerColor(kBlue);
        h2->SetMarkerStyle(21);
        h2->SetMarkerSize(0.8);
        h2->SetLineWidth(1);

        h3->SetLineColor(kRed);
        h3->SetMarkerColor(kRed);
        h3->SetMarkerStyle(22);
        h3->SetMarkerSize(0.8);
        h3->SetLineWidth(1);
        h3->Scale(10.0);  // Scale for visibility

        // Set titles and labels
        h1->SetTitle(titles[i]);
        h1->GetXaxis()->SetTitle("cos#theta_{GJ}");
        h1->GetYaxis()->SetTitle("Counts");
        h1->GetXaxis()->SetTitleSize(0.045);
        h1->GetYaxis()->SetTitleSize(0.045);
        h1->GetXaxis()->SetLabelSize(0.04);
        h1->GetYaxis()->SetLabelSize(0.04);

        // Find maximum for y-axis range
        double max1 = h1->GetMaximum();
        double max2 = h2->GetMaximum();
        double max3 = h3->GetMaximum();
        double ymax = std::max({max1, max2, max3}) * 1.2;

        h1->SetMaximum(ymax);

        // Draw histograms
        h1->Draw("HIST E1");
        h2->Draw("HIST E1 SAME");
        h3->Draw("HIST E1 SAME");

        // Add legend (only on first pad)
        if (i == 0) {
            TLegend* leg = new TLegend(0.18, 0.70, 0.45, 0.90);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(h1, "PPip", "lep");
            leg->AddEntry(h2, "Pip + P(FT)", "lep");
            leg->AddEntry(h3, "PPip + n cut (#times10)", "lep");
            leg->Draw();
        }
    }

    c1->Update();
    // c1->SaveAs("plot_pwa_groupE.pdf");
    // c1->SaveAs("plot_pwa_groupE.png");

    std::cout << "Plot displayed successfully" << std::endl;
}
