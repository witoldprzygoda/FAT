// isbest_fraction_betanew.C
//
// Data-only check of the isBest selection on the betanew production:
// for every pattern (e+e+, e-e-, e+e-) and trigger (PT3, PT2) the fraction
// of pairs with isBest==1 among all stored combinations ("isBest>=0";
// isBest==-1 rows are never written by FAT) as a function of the pair mass,
// and the like-sign CB = 2 sqrt(N++ N--) built from isBest==1 pairs versus
// from all pairs.  Exact identity used: CB_best / CB_all = sqrt(eps_pp eps_mm).
// The classic all / CB / signal plot is drawn for both variants (all stored
// combinations, i.e. isBest>=0, and isBest==1 with its own like-sign CB) and
// the two signals S = all - CB are compared directly.
// Also counts, per pair, the number of hypotheses of the same pattern in its
// event (adjacent rows with the same runnumber/event), because isBest==1 is
// automatic when the event has a single hypothesis.
//
// Input : output_epem_exp.root, output_epep_exp.root, output_emem_exp.root
//         (dilepton_nt, all branches Float_t: m_ee, trigbit, isBest, runnumber, event)
// Output: plots/output/isbest_fraction_betanew.root (all histograms),
//         plots/output/isbest_fraction_betanew_{PT3,PT2}.{pdf,png},
//         printed tables (per mass range and per 20 MeV MID bin).
//
// Errors: fractions carry plain binomial errors (pairs of one event are
// correlated through the ranking, so they are 20-45% optimistic for the
// like-sign fractions); CB_best/CB_all uses first-order propagation of
// sqrt(eps_pp eps_mm); signal errors add N+- and CB in quadrature.
//
// Run:  root -l -b -q 'plots/isbest_fraction_betanew.C+'

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TTree.h>
#include <TBranch.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

constexpr int    kNb   = 58;
constexpr double kMdLo = 0.14;
constexpr double kMdHi = 0.80;
constexpr float  kPT3  = 8192.f;
constexpr float  kPT2  = 4096.f;

const double* getEdges() {
    static double edges[kNb + 1];
    static bool built = false;
    if (!built) {
        int k = 0;
        for (int i = 0; i <= 20; ++i) edges[k++] = i * 0.0025;   // 0 .. 0.05
        edges[k++] = 0.0800; edges[k++] = 0.1100; edges[k++] = 0.1400;
        for (int i = 1; i <= 33; ++i) edges[k++] = 0.140 + i * 0.020;  // MID 20 MeV
        edges[k++] = 0.950; edges[k++] = 1.400;
        built = true;
    }
    return edges;
}

struct PatHists {
    // [trigger 0=PT3,1=PT2]
    TH1D* all[2];      // all stored pairs
    TH1D* best[2];     // isBest==1
    TH1D* one[2];      // pairs from events with exactly one hypothesis of this pattern
    TH1D* multi[2];    // pairs from events with >=2 hypotheses
    TH1D* multiBest[2];// isBest==1 among those
    TH1D* nhyp[2];     // hypotheses per event (one entry per event)
};

PatHists makeHists(const char* pat) {
    PatHists H;
    const char* trg[2] = {"PT3", "PT2"};
    for (int t = 0; t < 2; ++t) {
        H.all[t]       = new TH1D(Form("h_all_%s_%s", pat, trg[t]),       ";M_{ee} [GeV/c^{2}];pairs", kNb, getEdges());
        H.best[t]      = new TH1D(Form("h_best_%s_%s", pat, trg[t]),      ";M_{ee} [GeV/c^{2}];pairs", kNb, getEdges());
        H.one[t]       = new TH1D(Form("h_onehyp_%s_%s", pat, trg[t]),    ";M_{ee} [GeV/c^{2}];pairs", kNb, getEdges());
        H.multi[t]     = new TH1D(Form("h_multihyp_%s_%s", pat, trg[t]),  ";M_{ee} [GeV/c^{2}];pairs", kNb, getEdges());
        H.multiBest[t] = new TH1D(Form("h_multihyp_best_%s_%s", pat, trg[t]), ";M_{ee} [GeV/c^{2}];pairs", kNb, getEdges());
        H.nhyp[t]      = new TH1D(Form("h_nhyp_%s_%s", pat, trg[t]),      ";hypotheses per event;events", 12, 0.5, 12.5);
        for (TH1D* h : {H.all[t], H.best[t], H.one[t], H.multi[t], H.multiBest[t], H.nhyp[t]}) h->Sumw2();
    }
    return H;
}

// One pass over dilepton_nt. Rows of one event are adjacent (FAT writes the
// PAT hypotheses in order), so the hypothesis count is a run length.
Long64_t fillPattern(const char* file, PatHists& H) {
    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return -1; }
    TTree* t = (TTree*)f->Get("dilepton_nt");
    if (!t) { printf("no dilepton_nt in %s\n", file); return -1; }
    t->SetBranchStatus("*", 0);
    Float_t m = 0, trig = 0, best = 0, run = 0, ev = 0;
    for (const char* b : {"m_ee", "trigbit", "isBest", "runnumber", "event"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("m_ee", &m);
    t->SetBranchAddress("trigbit", &trig);
    t->SetBranchAddress("isBest", &best);
    t->SetBranchAddress("runnumber", &run);
    t->SetBranchAddress("event", &ev);

    const Long64_t N = t->GetEntries();
    std::vector<float> em, eb; std::vector<int> et;   // rows of the current event
    float curRun = -1.f, curEv = -1.f; bool have = false;
    Long64_t nev = 0, nbad = 0;

    auto flush = [&]() {
        if (em.empty()) return;
        const int nh = (int)em.size();
        int nbest = 0;
        for (size_t i = 0; i < em.size(); ++i) nbest += (eb[i] > 0.5f);
        if (nbest != 1) ++nbad;
        const int t0 = et[0];
        if (t0 >= 0) H.nhyp[t0]->Fill(nh);
        for (size_t i = 0; i < em.size(); ++i) {
            const int tt = et[i];
            if (tt < 0) continue;
            H.all[tt]->Fill(em[i]);
            if (eb[i] > 0.5f) H.best[tt]->Fill(em[i]);
            if (nh == 1) H.one[tt]->Fill(em[i]);
            else { H.multi[tt]->Fill(em[i]); if (eb[i] > 0.5f) H.multiBest[tt]->Fill(em[i]); }
        }
        ++nev;
        em.clear(); eb.clear(); et.clear();
    };

    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (!have || run != curRun || ev != curEv) { flush(); curRun = run; curEv = ev; have = true; }
        const int tt = (trig == kPT3) ? 0 : (trig == kPT2) ? 1 : -1;
        em.push_back(m); eb.push_back(best); et.push_back(tt);
        if ((i % 20000000) == 0 && i) printf("   %s: %lld / %lld rows\n", file, i, N);
    }
    flush();
    printf("%s: rows %lld, events %lld, events with #isBest!=1: %lld\n", file, N, nev, nbad);
    f->Close();
    return N;
}

double integ(TH1D* h, double lo, double hi) {
    const int b1 = h->FindBin(lo + 1e-6), b2 = h->FindBin(hi - 1e-6);
    return h->Integral(b1, b2);
}

// fraction best/all with binomial error
void frac(double nb, double na, double& e, double& de) {
    e = (na > 0) ? nb / na : 0.;
    de = (na > 0) ? std::sqrt(std::max(e * (1. - e), 1. / na) / na) : 0.;   // floor 1/N at eps = 0 or 1
}

TGraphErrors* fracGraph(TH1D* hb, TH1D* ha, int color, int marker) {
    TGraphErrors* g = new TGraphErrors();
    for (int b = 1; b <= ha->GetNbinsX(); ++b) {
        const double na = ha->GetBinContent(b), nb = hb->GetBinContent(b);
        if (na <= 0) continue;
        double e, de; frac(nb, na, e, de);
        const int n = g->GetN();
        g->SetPoint(n, ha->GetBinCenter(b), e);
        g->SetPointError(n, ha->GetBinWidth(b) / 2., de);
    }
    g->SetMarkerColor(color); g->SetLineColor(color); g->SetMarkerStyle(marker); g->SetMarkerSize(0.9);
    return g;
}

TGraphErrors* sqrtFracGraph(TH1D* hbp, TH1D* hap, TH1D* hbm, TH1D* ham, int color, int marker) {
    TGraphErrors* g = new TGraphErrors();
    for (int b = 1; b <= hap->GetNbinsX(); ++b) {
        const double nap = hap->GetBinContent(b), nbp = hbp->GetBinContent(b);
        const double nam = ham->GetBinContent(b), nbm = hbm->GetBinContent(b);
        if (nap <= 0 || nam <= 0) continue;
        double ep, dep, emm, dem; frac(nbp, nap, ep, dep); frac(nbm, nam, emm, dem);
        if (ep <= 0 || emm <= 0) continue;
        const double r = std::sqrt(ep * emm);
        const double dr = 0.5 * r * std::sqrt(std::pow(dep / ep, 2) + std::pow(dem / emm, 2));
        const int n = g->GetN();
        g->SetPoint(n, hap->GetBinCenter(b), r);
        g->SetPointError(n, hap->GetBinWidth(b) / 2., dr);
    }
    g->SetMarkerColor(color); g->SetLineColor(color); g->SetMarkerStyle(marker); g->SetMarkerSize(0.9);
    return g;
}

// CB_best / CB_all = sqrt(eps_pp eps_mm) with first-order error propagation (e+e+ and e-e- independent).
TH1D* sqrtFracHist(TH1D* hbp, TH1D* hap, TH1D* hbm, TH1D* ham, const char* name) {
    TH1D* r = (TH1D*)hap->Clone(name); r->SetDirectory(nullptr); r->Reset();
    for (int b = 1; b <= hap->GetNbinsX(); ++b) {
        const double nap = hap->GetBinContent(b), nbp = hbp->GetBinContent(b);
        const double nam = ham->GetBinContent(b), nbm = hbm->GetBinContent(b);
        if (nap <= 0 || nam <= 0) continue;
        double ep, dep, emm, dem; frac(nbp, nap, ep, dep); frac(nbm, nam, emm, dem);
        if (ep <= 0 || emm <= 0) continue;
        const double v = std::sqrt(ep * emm);
        r->SetBinContent(b, v);
        r->SetBinError(b, 0.5 * v * std::sqrt(std::pow(dep / ep, 2) + std::pow(dem / emm, 2)));
    }
    return r;
}

TH1D* makeCB(TH1D* hp, TH1D* hm, const char* name) {
    TH1D* cb = (TH1D*)hp->Clone(name); cb->SetDirectory(nullptr); cb->Reset();
    for (int b = 1; b <= hp->GetNbinsX(); ++b) {
        const double np = hp->GetBinContent(b), nm = hm->GetBinContent(b);
        if (np <= 0 || nm <= 0) { cb->SetBinContent(b, 0); cb->SetBinError(b, 0); continue; }
        const double v = 2. * std::sqrt(np * nm);
        cb->SetBinContent(b, v);
        cb->SetBinError(b, 0.5 * v * std::sqrt(1. / np + 1. / nm));
    }
    return cb;
}

// S = N+- - CB, errors added in quadrature (N+- and the like-sign counts are independent).
TH1D* makeSignal(TH1D* hn, TH1D* hcb, const char* name) {
    TH1D* s = (TH1D*)hn->Clone(name); s->SetDirectory(nullptr);
    for (int b = 1; b <= hn->GetNbinsX(); ++b) {
        s->SetBinContent(b, hn->GetBinContent(b) - hcb->GetBinContent(b));
        s->SetBinError(b, std::sqrt(std::pow(hn->GetBinError(b), 2) + std::pow(hcb->GetBinError(b), 2)));
    }
    return s;
}

void style(TH1D* h, int color, int marker) {
    h->SetMarkerColor(color); h->SetLineColor(color); h->SetMarkerStyle(marker); h->SetMarkerSize(0.9);
    h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.2);
}

} // namespace

void isbest_fraction_betanew(const char* file_pm = "output_epem_exp.root",
                             const char* file_pp = "output_epep_exp.root",
                             const char* file_mm = "output_emem_exp.root",
                             const char* out_dir = "plots/output")
{
    gStyle->SetOptStat(0);
    PatHists Hpp = makeHists("pp"), Hmm = makeHists("mm"), Hpm = makeHists("pm");
    if (fillPattern(file_pp, Hpp) < 0) return;
    if (fillPattern(file_mm, Hmm) < 0) return;
    if (fillPattern(file_pm, Hpm) < 0) return;

    const char* trg[2] = {"PT3", "PT2"};
    const double rlo[4] = {0.00, 0.14, 0.80, 0.00};
    const double rhi[4] = {0.14, 0.80, 1.40, 1.40};
    const char* rname[4] = {"M < 0.14", "0.14-0.80", "0.80-1.40", "0.00-1.40"};

    TFile* fout = TFile::Open(Form("%s/isbest_fraction_betanew.root", out_dir), "RECREATE");

    for (int t = 0; t < 2; ++t) {
        printf("\n================ %s ================\n", trg[t]);
        printf("hypotheses per event (events): ");
        for (const char* p : {"pp", "mm", "pm"}) {
            TH1D* hn = (p[0] == 'p' && p[1] == 'p') ? Hpp.nhyp[t] : (p[0] == 'm') ? Hmm.nhyp[t] : Hpm.nhyp[t];
            printf(" %s:", p);
            for (int b = 1; b <= 6; ++b) printf(" n=%d:%.0f", b, hn->GetBinContent(b));
            printf(" n>6:%.0f |", hn->Integral(7, 13));
        }
        printf("\n");

        printf("\n-- per mass range: pairs all / isBest==1 (fraction), share of pairs from 1-hypothesis events, fraction isBest==1 inside >=2-hypothesis events --\n");
        for (int r = 0; r < 4; ++r) {
            printf(" %-10s", rname[r]);
            for (const char* p : {"pp", "mm", "pm"}) {
                PatHists& H = (p[0] == 'p' && p[1] == 'p') ? Hpp : (p[0] == 'm') ? Hmm : Hpm;
                const double na = integ(H.all[t], rlo[r], rhi[r]), nb = integ(H.best[t], rlo[r], rhi[r]);
                const double n1 = integ(H.one[t], rlo[r], rhi[r]);
                const double nm = integ(H.multi[t], rlo[r], rhi[r]), nmb = integ(H.multiBest[t], rlo[r], rhi[r]);
                double e, de; frac(nb, na, e, de);
                printf(" | %s: %9.0f / %9.0f (%.4f+-%.4f)  1hyp %.3f  eps(>=2hyp) %.3f", p, na, nb, e, de, na > 0 ? n1 / na : 0., nm > 0 ? nmb / nm : 0.);
            }
            printf("\n");
        }
        printf("\n-- per mass range: CB = 2 sqrt(N++ N--) from all pairs vs isBest==1, ratio = sqrt(eps_pp eps_mm), compared with eps_pm; naive signals --\n");
        printf(" %-10s | %10s %10s %8s | %8s | %10s %10s | %10s %10s %8s\n", "range", "CB_all", "CB_best", "ratio", "eps_pm", "N+-_all", "N+-_best", "S_all", "S_best^nv", "S ratio");
        for (int r = 0; r < 4; ++r) {
            const double pa = integ(Hpp.all[t], rlo[r], rhi[r]), pb = integ(Hpp.best[t], rlo[r], rhi[r]);
            const double ma = integ(Hmm.all[t], rlo[r], rhi[r]), mb = integ(Hmm.best[t], rlo[r], rhi[r]);
            const double xa = integ(Hpm.all[t], rlo[r], rhi[r]), xb = integ(Hpm.best[t], rlo[r], rhi[r]);
            const double cba = 2. * std::sqrt(pa * ma), cbb = 2. * std::sqrt(pb * mb);
            const double sa = xa - cba, sb = xb - cbb;
            printf(" %-10s | %10.0f %10.0f %8.4f | %8.4f | %10.0f %10.0f | %10.0f %10.0f %8.3f\n",
                   rname[r], cba, cbb, cba > 0 ? cbb / cba : 0., xa > 0 ? xb / xa : 0., xa, xb, sa, sb, sa != 0 ? sb / sa : 0.);
        }

        printf("\n-- MID bins (20 MeV): eps_pp  eps_mm  sqrt(eps_pp eps_mm)  eps_pm | CB_all CB_best | N+-_all N+-_best | S_all S_best^naive  S_best/S_all --\n");
        for (int b = Hpp.all[t]->FindBin(kMdLo + 1e-6); b <= Hpp.all[t]->FindBin(kMdHi - 1e-6); ++b) {
            const double pa = Hpp.all[t]->GetBinContent(b), pb = Hpp.best[t]->GetBinContent(b);
            const double ma = Hmm.all[t]->GetBinContent(b), mb = Hmm.best[t]->GetBinContent(b);
            const double xa = Hpm.all[t]->GetBinContent(b), xb = Hpm.best[t]->GetBinContent(b);
            double ep, dep, em, dem, ex, dex; frac(pb, pa, ep, dep); frac(mb, ma, em, dem); frac(xb, xa, ex, dex);
            const double cba = 2. * std::sqrt(pa * ma), cbb = 2. * std::sqrt(pb * mb);
            printf("  %.2f-%.2f: %.3f+-%.3f  %.3f+-%.3f  %.3f  %.3f+-%.3f | %8.0f %8.0f | %8.0f %8.0f | %7.0f %7.0f  %.3f\n",
                   Hpp.all[t]->GetBinLowEdge(b), Hpp.all[t]->GetBinLowEdge(b + 1),
                   ep, dep, em, dem, std::sqrt(ep * em), ex, dex, cba, cbb, xa, xb, xa - cba, xb - cbb,
                   (xa - cba) != 0. ? (xb - cbb) / (xa - cba) : 0.);
        }

        // ---------- figure: 2 x 3 ----------
        TH1D* hxa = (TH1D*)Hpm.all[t]->Clone(Form("sp_pm_all_%s", trg[t]));  style(hxa, kBlack, 20);
        TH1D* hxb = (TH1D*)Hpm.best[t]->Clone(Form("sp_pm_best_%s", trg[t])); style(hxb, kBlack, 24);
        TH1D* cba = makeCB(Hpp.all[t], Hmm.all[t], Form("cb_all_%s", trg[t]));   style(cba, kRed + 1, 21);
        TH1D* cbb = makeCB(Hpp.best[t], Hmm.best[t], Form("cb_best_%s", trg[t])); style(cbb, kRed + 1, 25);
        TH1D* sga = makeSignal(hxa, cba, Form("sig_all_%s", trg[t]));  style(sga, kBlue + 1, 22);
        TH1D* sgb = makeSignal(hxb, cbb, Form("sig_best_%s", trg[t])); style(sgb, kBlue + 1, 26);
        for (TH1D* h : {hxa, hxb, cba, cbb, sga, sgb}) { h->SetTitle(";M_{ee} [GeV/c^{2}];pairs / bin"); h->GetXaxis()->SetRangeUser(0., 1.4); }
        const double ymax = 3. * hxa->GetMaximum();
        const double ymin = 0.5 * std::max(1., std::min(cba->GetMinimum(1e-9), sga->GetBinContent(sga->FindBin(0.79))));

        TCanvas* c = new TCanvas(Form("c_%s", trg[t]), "", 1800, 1100);
        c->Divide(3, 2);
        auto pad = [&](int i, bool logy) { c->cd(i); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.04); if (logy) gPad->SetLogy(); };
        TLatex lt; lt.SetNDC(); lt.SetTextSize(0.045);

        // (1) all stored combinations: all / CB / signal
        pad(1, true);
        hxa->SetMinimum(ymin); hxa->SetMaximum(ymax); hxa->Draw("E"); cba->Draw("E same"); sga->Draw("E same");
        TLegend* l1 = new TLegend(0.50, 0.60, 0.94, 0.82); l1->SetBorderSize(0); l1->SetTextSize(0.04);
        l1->AddEntry(hxa, "all: N_{+-}, all combinations", "p"); l1->AddEntry(cba, "CB = 2#sqrt{N_{++}N_{--}}, all comb.", "p"); l1->AddEntry(sga, "signal = all - CB", "p"); l1->Draw();
        lt.DrawLatex(0.16, 0.86, Form("%s, isBest #geq 0 (all stored combinations)", trg[t]));

        // (2) isBest==1: all / CB / signal
        pad(2, true);
        hxb->SetMinimum(ymin); hxb->SetMaximum(ymax); hxb->Draw("E"); cbb->Draw("E same"); sgb->Draw("E same");
        TLegend* l2 = new TLegend(0.50, 0.60, 0.94, 0.82); l2->SetBorderSize(0); l2->SetTextSize(0.04);
        l2->AddEntry(hxb, "all: N_{+-}, isBest = 1", "p"); l2->AddEntry(cbb, "CB = 2#sqrt{N_{++}N_{--}}, isBest = 1", "p"); l2->AddEntry(sgb, "signal = all - CB", "p"); l2->Draw();
        lt.DrawLatex(0.16, 0.86, Form("%s, isBest = 1 (one pair per event and pattern)", trg[t]));

        // (3) the two signals compared
        pad(3, true);
        TH1D* sga2 = (TH1D*)sga->Clone(Form("sig_all2_%s", trg[t])); sga2->SetMinimum(ymin); sga2->SetMaximum(ymax);
        sga2->Draw("E"); sgb->Draw("E same");
        TLegend* l3 = new TLegend(0.40, 0.64, 0.94, 0.82); l3->SetBorderSize(0); l3->SetTextSize(0.04);
        l3->AddEntry(sga2, "signal, all combinations", "p"); l3->AddEntry(sgb, "signal, isBest = 1 (naive CB)", "p"); l3->Draw();
        lt.DrawLatex(0.16, 0.86, Form("%s, signal = all - CB in both variants", trg[t]));

        // (4) fractions isBest==1 / all
        pad(4, false); gPad->SetGridy();
        TH1D* fr = new TH1D(Form("fr_%s", trg[t]), ";M_{ee} [GeV/c^{2}];fraction of pairs with isBest = 1", 10, 0., 1.4);
        fr->SetMinimum(0.5); fr->SetMaximum(1.05); style(fr, kBlack, 20); fr->Draw();
        TGraphErrors* gpp = fracGraph(Hpp.best[t], Hpp.all[t], kRed + 1, 20);
        TGraphErrors* gmm = fracGraph(Hmm.best[t], Hmm.all[t], kBlue + 1, 21);
        TGraphErrors* gpm = fracGraph(Hpm.best[t], Hpm.all[t], kBlack, 22);
        TGraphErrors* gsq = sqrtFracGraph(Hpp.best[t], Hpp.all[t], Hmm.best[t], Hmm.all[t], kGreen + 2, 24);
        gpp->Draw("P same"); gmm->Draw("P same"); gsq->Draw("P same"); gpm->Draw("P same");
        TLegend* l4 = new TLegend(0.16, 0.16, 0.62, 0.42); l4->SetBorderSize(0); l4->SetTextSize(0.04);
        l4->AddEntry(gpp, "e^{+}e^{+}: #varepsilon_{++}", "p"); l4->AddEntry(gmm, "e^{-}e^{-}: #varepsilon_{--}", "p");
        l4->AddEntry(gpm, "e^{+}e^{-}: #varepsilon_{+-}", "p"); l4->AddEntry(gsq, "#sqrt{#varepsilon_{++}#varepsilon_{--}} = CB_{best}/CB_{all}", "p"); l4->Draw();
        lt.DrawLatex(0.16, 0.86, Form("%s, isBest = 1 / all stored combinations", trg[t]));

        // (5) ratios isBest==1 / all for all, CB and signal
        pad(5, false); gPad->SetGridy();
        TH1D* fr5 = new TH1D(Form("fr5_%s", trg[t]), ";M_{ee} [GeV/c^{2}];isBest = 1 / all combinations", 10, 0., 1.4);
        fr5->SetMinimum(0.0); fr5->SetMaximum(1.1); style(fr5, kBlack, 20); fr5->Draw();
        TH1D* rpm = (TH1D*)hxb->Clone(Form("r_pm_%s", trg[t])); rpm->Divide(hxb, hxa, 1., 1., "B"); style(rpm, kBlack, 20);
        TH1D* rcb = sqrtFracHist(Hpp.best[t], Hpp.all[t], Hmm.best[t], Hmm.all[t], Form("r_cb_%s", trg[t])); style(rcb, kRed + 1, 21);
        TH1D* rsg = (TH1D*)sgb->Clone(Form("r_sig_%s", trg[t])); rsg->Divide(sgb, sga);            style(rsg, kBlue + 1, 22);
        rpm->Draw("E same"); rcb->Draw("E same"); rsg->Draw("E same");
        TLegend* l5 = new TLegend(0.16, 0.16, 0.62, 0.36); l5->SetBorderSize(0); l5->SetTextSize(0.04);
        l5->AddEntry(rpm, "all: N_{+-}^{best} / N_{+-}^{all}", "p"); l5->AddEntry(rcb, "CB_{best} / CB_{all}", "p"); l5->AddEntry(rsg, "signal_{best} / signal_{all}", "p"); l5->Draw();
        lt.DrawLatex(0.16, 0.86, Form("%s, ratios isBest = 1 / all combinations", trg[t]));

        // (6) like-sign spectra in both variants
        pad(6, true);
        TH1D* hpa = (TH1D*)Hpp.all[t]->Clone(Form("sp_pp_all_%s", trg[t]));  style(hpa, kRed + 1, 20);
        TH1D* hpb = (TH1D*)Hpp.best[t]->Clone(Form("sp_pp_best_%s", trg[t])); style(hpb, kRed + 1, 24);
        TH1D* hma = (TH1D*)Hmm.all[t]->Clone(Form("sp_mm_all_%s", trg[t]));  style(hma, kBlue + 1, 21);
        TH1D* hmb = (TH1D*)Hmm.best[t]->Clone(Form("sp_mm_best_%s", trg[t])); style(hmb, kBlue + 1, 25);
        for (TH1D* h : {hpa, hpb, hma, hmb}) { h->SetTitle(";M_{ee} [GeV/c^{2}];pairs / bin"); h->GetXaxis()->SetRangeUser(0., 1.4); }
        hma->SetMinimum(ymin); hma->SetMaximum(ymax); hma->Draw("E"); hmb->Draw("E same"); hpa->Draw("E same"); hpb->Draw("E same");
        TLegend* l6 = new TLegend(0.50, 0.56, 0.94, 0.82); l6->SetBorderSize(0); l6->SetTextSize(0.04);
        l6->AddEntry(hpa, "e^{+}e^{+}, all combinations", "p"); l6->AddEntry(hpb, "e^{+}e^{+}, isBest = 1", "p");
        l6->AddEntry(hma, "e^{-}e^{-}, all combinations", "p"); l6->AddEntry(hmb, "e^{-}e^{-}, isBest = 1", "p"); l6->Draw();
        lt.DrawLatex(0.16, 0.86, Form("%s, like-sign pairs", trg[t]));

        c->SaveAs(Form("%s/isbest_fraction_betanew_%s.pdf", out_dir, trg[t]));
        c->SaveAs(Form("%s/isbest_fraction_betanew_%s.png", out_dir, trg[t]));

        fout->cd();
        for (PatHists* H : {&Hpp, &Hmm, &Hpm})
            for (TH1D* h : {H->all[t], H->best[t], H->one[t], H->multi[t], H->multiBest[t], H->nhyp[t]}) h->Write();
        for (TH1D* h : {cba, cbb, sga, sgb, rcb, rpm, rsg}) h->Write();
        gpp->Write(Form("g_eps_pp_%s", trg[t])); gmm->Write(Form("g_eps_mm_%s", trg[t]));
        gpm->Write(Form("g_eps_pm_%s", trg[t])); gsq->Write(Form("g_sqrt_eps_%s", trg[t]));
    }
    fout->Close();
    printf("\nwritten: %s/isbest_fraction_betanew.root and isbest_fraction_betanew_{PT3,PT2}.{pdf,png}\n", out_dir);
}
