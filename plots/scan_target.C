// scan_target.C -- vertex picture of the target region from data.
// One row per event (isBest == 1) of the pi+pi-e+e- sample; like-sign samples for the combinatorial background.
// Input: the compact scan file written by scan_prepare.C. Output: plots/output/eta_scan/target_vertex_view.png
// Run: root -l -b -q plots/scan_target.C
void scan_target() {
    gStyle->SetOptStat(0); gStyle->SetPalette(kRainBow); gStyle->SetNumberContours(99);
    TFile f("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan/scan_compact.root");
    TTree* us = (TTree*)f.Get("us"); TTree* pp = (TTree*)f.Get("pp"); TTree* mm = (TTree*)f.Get("mm");
    const char* ev = "isBest>0.5"; const char* rt = "sqrt(vx*vx+vy*vy)";
    TH1D* hz = new TH1D("hz", "", 260, -260, 0); TH1D* hp = new TH1D("hp", "", 52, -260, 0); TH1D* hm = new TH1D("hm", "", 52, -260, 0); TH1D* hz5 = new TH1D("hz5", "", 52, -260, 0);
    us->Draw("vz>>hz", ev, "goff"); us->Draw("vz>>hz5", ev, "goff"); pp->Draw("vz>>hp", ev, "goff"); mm->Draw("vz>>hm", ev, "goff");
    double plat = hz->GetMaximum(), zlo = 0, zhi = 0; const int bmax = hz->GetMaximumBin();
    for (int b = bmax; b >= 2; --b) if (hz->GetBinContent(b - 1) < 0.5 * plat) { const double y1 = hz->GetBinContent(b - 1), y2 = hz->GetBinContent(b); zlo = hz->GetBinCenter(b - 1) + (0.5 * plat - y1) / (y2 - y1); break; }
    for (int b = bmax; b < 260; ++b) if (hz->GetBinContent(b + 1) < 0.5 * plat) { const double y1 = hz->GetBinContent(b), y2 = hz->GetBinContent(b + 1); zhi = hz->GetBinCenter(b) + (y1 - 0.5 * plat) / (y1 - y2); break; }
    printf("half-maximum edges: %.1f and %.1f mm, length %.1f mm; maximum at %.1f mm\n", zlo, zhi, zhi - zlo, hz->GetBinCenter(bmax));
    TH1D* hcb = (TH1D*)hz5->Clone("hcb"); hcb->Reset(); TH1D* hsh = (TH1D*)hz5->Clone("hsh"); hsh->Reset(); TH1D* hsig = (TH1D*)hz5->Clone("hsig"); hsig->Reset();
    for (int b = 1; b <= 52; ++b) { const double a = hp->GetBinContent(b), c = hm->GetBinContent(b), n = hz5->GetBinContent(b), cb = 2 * sqrt(a * c), ecb = sqrt(a + c);
        hcb->SetBinContent(b, cb); hcb->SetBinError(b, ecb); hsig->SetBinContent(b, n - cb); hsig->SetBinError(b, sqrt(n + a + c)); if (n > 0) { hsh->SetBinContent(b, cb / n); hsh->SetBinError(b, ecb / n); } }
    TH2D* hzr = new TH2D("hzr", "", 130, -260, 0, 60, 0, 30); us->Draw(Form("%s:vz>>hzr", rt), ev, "goff");
    TH2D* hxy = new TH2D("hxy", "", 80, -20, 20, 80, -20, 20); us->Draw("vy:vx>>hxy", Form("%s && vz>-160 && vz<-95", ev), "goff");
    const int NR = 4; const double rlo[NR] = {0, 8, 12, 18}, rhi[NR] = {4, 12, 18, 30}; const int col[NR] = {kBlack, kBlue, kGreen + 2, kRed}; TH1D* hr[NR];
    for (int i = 0; i < NR; ++i) { hr[i] = new TH1D(Form("hr%d", i), "", 80, -240, -40); us->Draw(Form("vz>>hr%d", i), Form("%s && %s>=%g && %s<%g", ev, rt, rlo[i], rt, rhi[i]), "goff"); hr[i]->Scale(1. / hr[i]->GetMaximum()); }
    TCanvas c("c", "", 2400, 1500); c.Divide(3, 2, 0.004, 0.004);
    auto style = [](TH1* h, int col, int ms, const char* xt, const char* yt) { h->SetMarkerStyle(ms); h->SetMarkerSize(0.9); h->SetMarkerColor(col); h->SetLineColor(col); h->SetLineWidth(2); h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
        h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05); h->GetXaxis()->SetLabelSize(0.042); h->GetYaxis()->SetLabelSize(0.042); h->GetYaxis()->SetTitleOffset(1.75); h->GetXaxis()->SetTitleOffset(1.1); h->GetXaxis()->SetNdivisions(506); };
    auto pad = [&](int i, bool logy = false, bool logz = false) { c.cd(i); gPad->SetLeftMargin(0.18); gPad->SetBottomMargin(0.14); gPad->SetRightMargin(logz ? 0.15 : 0.04); gPad->SetTopMargin(0.10); gPad->SetLogy(logy); gPad->SetLogz(logz); gPad->SetGrid(1, 1); };
    TLatex t; t.SetNDC(); t.SetTextSize(0.045); TLine l; l.SetLineWidth(2);
    pad(1); style(hz, kBlack, 20, "vertex z [mm]", "events / 1 mm"); hz->GetXaxis()->SetRangeUser(-240, -40); hz->SetMaximum(1.3 * plat); hz->Draw("E");
    l.SetLineStyle(2); l.SetLineColor(kRed + 1); l.DrawLine(-160, 0, -160, 1.2 * plat); l.SetLineColor(kGray + 2); l.SetLineStyle(3); l.DrawLine(zlo, 0, zlo, 1.05 * plat); l.DrawLine(zhi, 0, zhi, 1.05 * plat);
    t.DrawLatex(0.18, 0.93, "#pi^{+}#pi^{-}e^{+}e^{-} events, linear scale"); t.SetTextColor(kRed + 1); t.DrawLatex(0.21, 0.84, "z = -160"); t.SetTextColor(kGray + 3); t.DrawLatex(0.46, 0.84, Form("half maximum: %.0f, %.0f mm", zlo, zhi)); t.SetTextColor(kBlack);
    pad(2, true); style(hz5, kBlack, 20, "vertex z [mm]", "events / 5 mm"); hz5->SetMinimum(30); hz5->SetMaximum(3e6); hz5->Draw("E"); style(hcb, kRed, 21, "", ""); hcb->Draw("E same"); style(hsig, kBlue, 24, "", ""); hsig->Draw("E same");
    l.SetLineStyle(2); l.SetLineColor(kRed + 1); l.DrawLine(-160, 30, -160, 6e5); t.DrawLatex(0.18, 0.93, "logarithmic scale");
    TLegend lg(0.20, 0.72, 0.56, 0.89); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.042); lg.AddEntry(hz5, "all pairs", "p"); lg.AddEntry(hcb, "CB = 2#sqrt{N_{++}N_{--}}", "p"); lg.AddEntry(hsig, "all - CB", "p"); lg.Draw();
    pad(3); style(hsh, kRed, 21, "vertex z [mm]", "CB / all pairs"); hsh->GetXaxis()->SetRangeUser(-240, -40); hsh->SetMinimum(0); hsh->SetMaximum(0.6); hsh->Draw("E"); l.DrawLine(-160, 0, -160, 0.6); t.DrawLatex(0.18, 0.93, "share of CB (no k factor)");
    pad(4, false, true); style(hzr, kBlack, 20, "vertex z [mm]", "distance from the beam axis [mm]"); hzr->GetYaxis()->SetTitleOffset(1.3); hzr->Draw("COLZ"); l.SetLineColor(kBlack); l.DrawLine(-160, 0, -160, 30); t.DrawLatex(0.18, 0.93, "z against distance from the axis");
    pad(5); style(hr[0], col[0], 20, "vertex z [mm]", "events, each curve scaled to maximum 1"); hr[0]->SetMaximum(1.35); hr[0]->Draw("HIST"); for (int i = 1; i < NR; ++i) { style(hr[i], col[i], 20, "", ""); hr[i]->Draw("HIST same"); }
    l.SetLineColor(kRed + 1); l.DrawLine(-160, 0, -160, 1.1); t.DrawLatex(0.18, 0.93, "z in four ranges of that distance");
    TLegend lr(0.56, 0.68, 0.95, 0.89); lr.SetBorderSize(0); lr.SetFillStyle(0); lr.SetTextSize(0.042); for (int i = 0; i < NR; ++i) lr.AddEntry(hr[i], Form("%g - %g mm", rlo[i], rhi[i]), "l"); lr.Draw();
    pad(6, false, true); style(hxy, kBlack, 20, "vertex x [mm]", "vertex y [mm]"); hxy->GetYaxis()->SetTitleOffset(1.3); hxy->Draw("COLZ"); t.DrawLatex(0.18, 0.93, "x - y for -160 < z < -95 mm");
    c.SaveAs("/home/przygoda/HADES/FAT/FAT/plots/output/eta_scan/target_vertex_view.png");
}
