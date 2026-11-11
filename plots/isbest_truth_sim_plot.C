// isbest_truth_sim_plot.C
//
// Figure for the isBest-vs-MC-truth study (input: isbest_truth_sim_<tag>.root
// from plots/isbest_truth_sim.C). For events with >= 3 identified leptons
// (classes n = 3, 4, >= 5 summed; class 2 has no choice) it draws, versus the
// e+e- invariant mass:
//   left : eps_S = fraction of GENUINE pairs chosen
//   right: eps_B = fraction of COMBINATORIAL pairs chosen
// for the production ranking (tracking chi2, sim isBest), the time-chi2 ranking
// (the exp betanew isBest) and a random choice. If a ranking recognised the
// genuine pair, eps_S would lie above the random curve and eps_B below it.
//
// Run: root -l -b -q 'plots/isbest_truth_sim_plot.C+("std_pt3,lepAB_pt3","SMASH standard,SMASH lepton")'

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TStyle.h>
#include <vector>

namespace itsp {
// sum classes 3..5 and rebin to a coarse common grid for readable ratios
TH1D* sumClasses(TFile& f, const char* rank, const char* kind, const char* name) {
    TH1D* out = nullptr;
    for (int c = 3; c <= 5; ++c) {
        TH1D* h = (TH1D*)f.Get(Form("h_%s_%s_c%d", rank, kind, c));
        if (!h) continue;
        if (!out) { out = (TH1D*)h->Clone(name); out->SetDirectory(nullptr); out->Sumw2(); }
        else out->Add(h);
    }
    return out;
}
TH1D* rebinned(TH1D* h, const char* name) {
    const double e[] = {0.0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.11, 0.14, 0.20, 0.30, 0.40, 0.50, 0.65, 0.80, 1.40};
    const int n = sizeof(e) / sizeof(e[0]) - 1;
    TH1D* r = (TH1D*)h->Rebin(n, name, e); r->SetDirectory(nullptr);
    return r;
}
}

void isbest_truth_sim_plot(const char* tags = "std_pt3", const char* labels = "SMASH standard",
                           const char* in_dir = "plots/output", const char* out_name = "plots/output/isbest_truth_sim")
{
    using namespace itsp;
    gStyle->SetOptStat(0);
    TObjArray* T = TString(tags).Tokenize(",");
    TObjArray* L = TString(labels).Tokenize(",");
    const int nS = T->GetEntries();
    TCanvas* c = new TCanvas("c", "", 1500, 560 * nS);
    c->Divide(2, nS);
    const char* rk[3] = {"prod", "time", "rand"};
    const char* rkLab[3] = {"tracking #chi^{2} (sim isBest)", "time #chi^{2} (exp isBest)", "random choice"};
    const int col[3] = {kBlue + 1, kRed + 1, kGray + 2};
    const int mk[3] = {20, 21, 24};
    for (int s = 0; s < nS; ++s) {
        const TString tag = ((TObjString*)T->At(s))->GetString();
        const TString lab = s < L->GetEntries() ? ((TObjString*)L->At(s))->GetString() : tag;
        TFile* f = TFile::Open(Form("%s/isbest_truth_sim_%s.root", in_dir, tag.Data()));
        if (!f || f->IsZombie()) { printf("cannot open %s\n", tag.Data()); continue; }
        for (int k = 0; k < 2; ++k) {   // 0 genuine, 1 combinatorial
            const char* kind = k == 0 ? "gen" : "comb";
            c->cd(2 * s + k + 1);
            gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetGridy();
            TH1D* hA = rebinned(sumClasses(*f, "all", kind, Form("a_%d_%d", s, k)), Form("ar_%d_%d", s, k));
            TH1D* fr = new TH1D(Form("fr_%d_%d", s, k), "", 10, 0., 1.0);
            fr->SetMinimum(0.0); fr->SetMaximum(1.0);
            fr->GetXaxis()->SetTitle("M_{ee} [GeV/c^{2}]");
            fr->GetYaxis()->SetTitle(k == 0 ? "#varepsilon_{S}: genuine pairs chosen" : "#varepsilon_{B}: combinatorial pairs chosen");
            for (TAxis* a : {fr->GetXaxis(), fr->GetYaxis()}) { a->SetTitleSize(0.055); a->SetLabelSize(0.05); }
            fr->GetYaxis()->SetTitleOffset(1.1);
            fr->Draw();
            TLegend* lg = new TLegend(0.36, 0.17, 0.97, 0.37); lg->SetBorderSize(0); lg->SetTextSize(0.045); lg->SetFillStyle(0);
            for (int r = 0; r < 3; ++r) {
                TH1D* hS = rebinned(sumClasses(*f, rk[r], kind, Form("s_%d_%d_%d", s, k, r)), Form("sr_%d_%d_%d", s, k, r));
                TH1D* e = (TH1D*)hS->Clone(Form("e_%d_%d_%d", s, k, r));
                e->Divide(hS, hA, 1., 1., "B");
                for (int b = 1; b <= e->GetNbinsX(); ++b) if (hA->GetBinContent(b) < 30) { e->SetBinContent(b, -1); e->SetBinError(b, 0); }
                e->SetMarkerStyle(mk[r]); e->SetMarkerColor(col[r]); e->SetLineColor(col[r]); e->SetMarkerSize(1.1);
                e->Draw("E1 same");
                lg->AddEntry(e, rkLab[r], "p");
            }
            lg->Draw();
            TLatex t; t.SetNDC(); t.SetTextSize(0.05);
            t.DrawLatex(0.13, 0.92, Form("%s, PT3, events with #geq 3 leptons", lab.Data()));
        }
    }
    c->SaveAs(Form("%s.pdf", out_name));
    c->SaveAs(Form("%s.png", out_name));
}
