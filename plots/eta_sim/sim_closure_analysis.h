// ============================================================================
// sim_closure_analysis.h  --  stage 2 of sim_closure.C: tables and figures
// ============================================================================
// Included by sim_closure.C only (needs its namespace scl: binnings, names,
// constants). Reads the six pass files closure_pass_<k>.root, adds std and lep
// (independent event groups), writes
//   tables: plots/eta_sim/tables/closure_{ratios,ratio_bins,partition,flags,ghost,codes}.txt
//   figures: plots/eta_sim/figs/closure_{summary,pi0D,etaD,conv,modeA,partition_pre,partition_oa2,flags,ghost}.{png,pdf}
// ============================================================================
#ifndef SIM_CLOSURE_ANALYSIS_H
#define SIM_CLOSURE_ANALYSIS_H

namespace scl {

TFile* gF[6] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
const char* const kCGShort[kNCG] = {"#pi^{0} Dalitz pair", "#eta Dalitz pair", "two Dalitz decays", "conv. of #pi^{0} Dalitz #gamma",
                                    "conv. of #eta Dalitz #gamma", "mode A: std dropped vs lep owned"};

// ---- loaded accumulators
struct GV {
    std::vector<double> s, v, n;
    bool empty() const { return s.empty(); }
    double S(int i) const { return s.empty() ? 0. : s[i]; }
    double V(int i) const { return v.empty() ? 0. : v[i]; }
    double N(int i) const { return n.empty() ? 0. : n[i]; }
    double Neff(int i) const { return V(i) > 0 ? S(i) * S(i) / V(i) : 0.; }
    double Tot() const { return s.empty() ? 0. : s.back(); }
    double TotV() const { return v.empty() ? 0. : v.back(); }
    double TotN() const { return n.empty() ? 0. : n.back(); }
};
struct PV { std::vector<double> a, b, aa, bb, ab, n; bool empty() const { return a.empty(); } };

std::vector<double> Vec(TFile* f, const std::string& nm) {
    TVectorD* x = f->Get<TVectorD>(nm.c_str());
    if (!x) return {};
    std::vector<double> r(x->GetMatrixArray(), x->GetMatrixArray() + x->GetNrows());
    delete x;
    return r;
}
inline void AddVec(std::vector<double>& a, const std::vector<double>& b) {
    if (b.empty()) return;
    if (a.empty()) { a = b; return; }
    for (size_t i = 0; i < a.size(); ++i) a[i] += b[i];
}
GV LoadG(const std::string& base, int st, int h) {
    GV r;
    TFile* f = gF[3 * st + h];
    const std::string nm = base + Form("_t%d_h%d", st, h);
    r.s = Vec(f, nm + "__s"); r.v = Vec(f, nm + "__v"); r.n = Vec(f, nm + "__n");
    return r;
}
GV SumG(const std::string& base, const std::vector<std::pair<int, int>>& trees) {
    GV r;
    for (auto& t : trees) {
        GV x = LoadG(base, t.first, t.second);
        AddVec(r.s, x.s); AddVec(r.v, x.v); AddVec(r.n, x.n);
    }
    return r;
}
PV LoadP(const std::string& base, const std::vector<std::pair<int, int>>& trees) {
    PV r;
    for (auto& t : trees) {
        TFile* f = gF[3 * t.first + t.second];
        const std::string nm = base + Form("_t%d_h%d", t.first, t.second);
        AddVec(r.a, Vec(f, nm + "__a")); AddVec(r.b, Vec(f, nm + "__b")); AddVec(r.aa, Vec(f, nm + "__aa"));
        AddVec(r.bb, Vec(f, nm + "__bb")); AddVec(r.ab, Vec(f, nm + "__ab")); AddVec(r.n, Vec(f, nm + "__n"));
    }
    return r;
}
const std::vector<std::pair<int, int>> kStdLepUS = {{0, 0}, {1, 0}};

// ---- ratios
struct Rat { double r = 0, e = 0; };
inline Rat RatioIndep(double num, double vnum, double den, double vden) {
    Rat x;
    if (den > 0 && num > 0) { x.r = num / den; x.e = x.r * std::sqrt(vnum / (num * num) + vden / (den * den)); }
    return x;
}
// ratio estimator B/A over the same rows, group-correlated: V = sum_g (b_g - R a_g)^2 / A^2
inline Rat RatioPaired(const PV& p, int i) {
    Rat x;
    if (p.empty() || !(p.a[i] > 0)) return x;
    x.r = p.b[i] / p.a[i];
    const double V = p.bb[i] - 2 * x.r * p.ab[i] + x.r * x.r * p.aa[i];
    x.e = std::sqrt(std::max(0., V)) / p.a[i];
    return x;
}
std::string PM(double v, double e, const char* fmt = "%.4f") {
    std::string f = std::string(fmt) + " +- " + fmt;
    return Form(f.c_str(), v, e);
}
std::string Sig(double v, double e) {   // cross section with an adapted number of digits
    if (v == 0 && e == 0) return "0";
    const double m = std::max(std::fabs(v), e);
    const char* f = m >= 100 ? "%.1f" : m >= 10 ? "%.2f" : m >= 1 ? "%.3f" : m >= 0.1 ? "%.4f" : m >= 0.01 ? "%.5f" : "%.2e";
    return PM(v, e, f);
}

// ---- histograms
TH1D* HistOf(const GV& g, const Bins& b, const char* name, double scale = 1.) {   // dsigma/dx [nb / unit]
    TH1D* h = new TH1D(name, "", b.n(), b.e.data());
    h->SetDirectory(nullptr);
    h->Sumw2();
    for (int i = 0; i < b.n(); ++i) {
        const double wd = b.hi(i) - b.lo(i);
        const double f = scale / (simeng::kLsim * wd);
        h->SetBinContent(i + 1, g.S(i) * f);
        h->SetBinError(i + 1, std::sqrt(g.V(i)) * f);
    }
    return h;
}
TH1D* RatioHist(const GV& num, const GV& den, const Bins& b, const char* name, double neffMin) {
    TH1D* h = new TH1D(name, "", b.n(), b.e.data());
    h->SetDirectory(nullptr);
    for (int i = 0; i < b.n(); ++i) {
        if (num.Neff(i) < neffMin || den.Neff(i) < neffMin) { h->SetBinContent(i + 1, -999); h->SetBinError(i + 1, 0); continue; }
        const Rat r = RatioIndep(num.S(i), num.V(i), den.S(i), den.V(i));
        h->SetBinContent(i + 1, r.r);
        h->SetBinError(i + 1, r.e);
    }
    return h;
}
TH1D* RatioHistP(const PV& p, const Bins& b, const char* name) {
    TH1D* h = new TH1D(name, "", b.n(), b.e.data());
    h->SetDirectory(nullptr);
    for (int i = 0; i < b.n(); ++i) {
        const Rat r = RatioPaired(p, i);
        if (r.r == 0) { h->SetBinContent(i + 1, -999); continue; }
        h->SetBinContent(i + 1, r.r);
        h->SetBinError(i + 1, r.e);
    }
    return h;
}
// shape chi2 of two spectra (bins with N_eff >= neffMin in both), the second scaled to the first
struct Chi { double chi2 = 0; int ndf = 0; };
Chi ShapeChi2(const GV& a, const GV& b, int nb, double neffMin) {
    double sa = 0, sb = 0;
    std::vector<int> use;
    for (int i = 0; i < nb; ++i)
        if (a.Neff(i) >= neffMin && b.Neff(i) >= neffMin) { use.push_back(i); sa += a.S(i); sb += b.S(i); }
    Chi c;
    if (use.size() < 2 || sb <= 0) return c;
    const double k = sa / sb;
    for (int i : use) c.chi2 += std::pow(a.S(i) - k * b.S(i), 2) / (a.V(i) + k * k * b.V(i));
    c.ndf = (int)use.size() - 1;
    return c;
}

// ---- plotting helpers
void StyleLine(TH1D* h, int col, int lw = 2, int ls = 1, int fill = 0) {
    h->SetLineColor(col); h->SetLineWidth(lw); h->SetLineStyle(ls); h->SetMarkerSize(0); h->SetMarkerColor(col);
    if (fill) { h->SetFillColor(fill); h->SetFillStyle(1001); } else h->SetFillStyle(0);
}
void StyleMark(TH1D* h, int col, int mst, double ms = 1.1) {
    h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(mst); h->SetMarkerSize(ms); h->SetLineWidth(2);
}
void AxisStyle(TH1* h, const char* xt, const char* yt, double sc = 1.0) {
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(0.058 * sc); h->GetYaxis()->SetTitleSize(0.058 * sc);
    h->GetXaxis()->SetLabelSize(0.050 * sc); h->GetYaxis()->SetLabelSize(0.050 * sc);
    h->GetXaxis()->SetTitleOffset(1.05); h->GetYaxis()->SetTitleOffset(1.2 / sc);
    h->GetYaxis()->SetMaxDigits(4);
}
void Txt(double x, double y, const char* s, double size = 0.04, int col = 1, int align = 11) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextColor(col); t.SetTextAlign(align); t.DrawLatex(x, y, s);
}
void PadMargins(TVirtualPad* p, double l = 0.15, double r = 0.03, double t = 0.08, double b = 0.14) {
    p->SetLeftMargin(l); p->SetRightMargin(r); p->SetTopMargin(t); p->SetBottomMargin(b);
}
void Save(TCanvas* c, const std::string& base, bool pdf) {
    c->SaveAs((std::string(simeng::kFigDir) + "/" + base + ".png").c_str());
    if (pdf) c->SaveAs((std::string(simeng::kFigDir) + "/" + base + ".pdf").c_str());
}
double MaxOf(const std::vector<TH1D*>& hs, double xlo, double xhi) {
    double m = 0;
    for (TH1D* h : hs)
        for (int i = 1; i <= h->GetNbinsX(); ++i)
            if (h->GetBinCenter(i) > xlo && h->GetBinCenter(i) < xhi) m = std::max(m, h->GetBinContent(i) + h->GetBinError(i));
    return m;
}
double MinPos(const std::vector<TH1D*>& hs, double xlo, double xhi) {
    double m = 1e30;
    for (TH1D* h : hs)
        for (int i = 1; i <= h->GetNbinsX(); ++i)
            if (h->GetBinCenter(i) > xlo && h->GetBinCenter(i) < xhi && h->GetBinContent(i) > 0) m = std::min(m, h->GetBinContent(i));
    return m < 1e30 ? m : 1e-6;
}
void WindowBoxes(double ylo, double yhi, bool log) {
    // bump and peak windows (shaded), above-peak window (hatched outline)
    const double wl[3] = {kWinEdge[1], kWinEdge[3], kWinEdge[5]}, wh[3] = {kWinEdge[2], kWinEdge[4], kWinEdge[6]};
    const int col[3] = {kOrange - 9, kGreen - 9, kAzure - 9};
    for (int i = 0; i < 3; ++i) {
        TBox* b = new TBox(wl[i], ylo, wh[i], yhi);
        b->SetFillColorAlpha(col[i], 0.35); b->SetLineWidth(0);
        b->Draw();
    }
    (void)log;
}

// ============================================================================
// Part 1: closure
// ============================================================================
struct Audit { double r, e; const char* src; };
void PartClosure(FILE* ft, FILE* fb) {
    const char* hypN[3] = {"e+e-", "e+e+", "e-e-"};
    const char* selN[2] = {"pre", "pre+OA>2"};
    fprintf(ft, "# closure_ratios: lep/std of the rows the lepton set provides in mode A (std + lep FAT outputs, every row)\n");
    fprintf(ft, "# selection pre = isBest>=0, PT3, vz>-160, cor_ok; pre+OA>2 adds oa > 2 deg. Weights: sim_w (PDG BR) on both sides,\n");
    fprintf(ft, "#   conversions of Dalitz photons: std = every photon of a primary pi0 (eta) x BR_D(PDG)/N_gamma(HGeant) = %.6f (%.6f),\n", kFconvPi0, kFconvEta);
    fprintf(ft, "#   lep = Dalitz-photon conversion x BR_D(PDG) (0.01174, 0.0069). Expected ratio 1 everywhere.\n");
    fprintf(ft, "# pions from the common SMASH final state: no e+- and no rho0/omega/phi/eta product or forced-decay product in a pion slot\n");
    fprintf(ft, "# sigma = sum w / 221.05 nb^-1; errors: event-group variance, std and lep independent (shared SMASH events neglected)\n\n");
    fprintf(ft, "%-44s %-6s %-9s %8s %18s %9s %9s %18s %9s %18s %18s\n", "group", "pair", "sel", "std rows", "std sigma [nb]",
            "std Neff", "lep rows", "lep sigma [nb]", "lep Neff", "lep/std", "lep/std no ghost");
    for (int g = 0; g < kNCG; ++g)
        for (int h = 0; h < 3; ++h) {
            if (h > 0 && g != kDD && g != kModeA) continue;
            for (int s = 0; s < 2; ++s) {
                const GV a = LoadG(Form("cl_%s_s%d_o%d", kCGKey[g], s, kOTot), 0, h);
                const GV b = LoadG(Form("cl_%s_s%d_o%d", kCGKey[g], s, kOTot), 1, h);
                const GV a0 = LoadG(Form("cl_%s_s%d_o%d", kCGKey[g], s, kOTotNoGhost), 0, h);
                const GV b0 = LoadG(Form("cl_%s_s%d_o%d", kCGKey[g], s, kOTotNoGhost), 1, h);
                const Rat r = RatioIndep(b.Tot(), b.TotV(), a.Tot(), a.TotV());
                const Rat r0 = RatioIndep(b0.Tot(), b0.TotV(), a0.Tot(), a0.TotV());
                fprintf(ft, "%-44s %-6s %-9s %8.0f %18s %9.0f %9.0f %18s %9.0f %18s %18s\n", kCGName[g], hypN[h], selN[s], a.TotN(),
                        Sig(a.Tot() / simeng::kLsim, std::sqrt(a.TotV()) / simeng::kLsim).c_str(), a.Neff(a.s.size() - 1),
                        b.TotN(), Sig(b.Tot() / simeng::kLsim, std::sqrt(b.TotV()) / simeng::kLsim).c_str(),
                        b.Neff(b.s.size() - 1), PM(r.r, r.e).c_str(), PM(r0.r, r0.e).c_str());
            }
        }
    // two Dalitz decays, like-sign trees together
    for (int s = 0; s < 2; ++s) {
        const GV a = SumG(Form("cl_DD_s%d_o%d", s, kOTot), {{0, 1}, {0, 2}});
        const GV b = SumG(Form("cl_DD_s%d_o%d", s, kOTot), {{1, 1}, {1, 2}});
        const Rat r = RatioIndep(b.Tot(), b.TotV(), a.Tot(), a.TotV());
        fprintf(ft, "%-44s %-6s %-9s %8.0f %18s %9.0f %9.0f %18s %9.0f %18s\n", "two Dalitz decays", "ls", selN[s], a.TotN(),
                Sig(a.Tot() / simeng::kLsim, std::sqrt(a.TotV()) / simeng::kLsim).c_str(), a.Neff(a.s.size() - 1), b.TotN(),
                Sig(b.Tot() / simeng::kLsim, std::sqrt(b.TotV()) / simeng::kLsim).c_str(), b.Neff(b.s.size() - 1),
                PM(r.r, r.e).c_str());
    }
    fprintf(ft, "\n# audit values (audit_closure.md, PAT trees, PT3 + z>-500, no vz>-160 and no OA cut; sim_weights_spec.md closure table):\n");
    fprintf(ft, "#   pi0 Dalitz, common pions, PDG 0.894 +- 0.003 | eta Dalitz 0.931 +- 0.005 | two Dalitz e+e- 1.19 +- 0.18 (spec 1.04 +- 0.11),\n");
    fprintf(ft, "#   like-sign 1.09 +- 0.14 (spec e+e+ 0.94 +- 0.16, e-e- 1.00 +- 0.12) | conv. pi0 photon 0.952 +- 0.009 | conv. eta photon 1.06 +- 0.03\n");
    fprintf(ft, "#   mode A overall (rows lep provides / rows std drops): e+e- 0.896 +- 0.004, e+e+ 1.10 +- 0.14, e-e- 0.96 +- 0.08\n");

    // ---- per-bin ratios and shape chi2
    const char* obsN[4] = {"", "M(ee) [GeV/c2]", "OA [deg]", "M4 [GeV/c2]"};
    const Bins* obsB[4] = {&bTot, &bMee, &bOA, &bM4c};
    fprintf(fb, "# closure_ratio_bins: lep/std per bin (e+e-), bins with N_eff >= 10 in both sets; shape chi2 = std vs lep scaled to std\n");
    for (int g = 0; g < kNCG; ++g) {
        if (g == kDD) continue;
        for (int o = 1; o <= 3; ++o)
            for (int s = 0; s < 2; ++s) {
                const GV a = LoadG(Form("cl_%s_s%d_o%d", kCGKey[g], s, o), 0, 0);
                const GV b = LoadG(Form("cl_%s_s%d_o%d", kCGKey[g], s, o), 1, 0);
                const Chi c = ShapeChi2(a, b, obsB[o]->n(), 10.);
                fprintf(fb, "\n## %s, %s, %s: shape chi2/ndf = %.1f/%d\n", kCGName[g], selN[s], obsN[o], c.chi2, c.ndf);
                fprintf(fb, "%10s %10s %10s %10s %10s %18s\n", "lo", "hi", "std Neff", "lep Neff", "std rows", "lep/std");
                for (int i = 0; i < obsB[o]->n(); ++i) {
                    if (a.Neff(i) < 10 || b.Neff(i) < 10) continue;
                    const Rat r = RatioIndep(b.S(i), b.V(i), a.S(i), a.V(i));
                    fprintf(fb, "%10.4g %10.4g %10.0f %10.0f %10.0f %18s\n", obsB[o]->lo(i), obsB[o]->hi(i), a.Neff(i), b.Neff(i),
                            a.N(i), PM(r.r, r.e).c_str());
                }
            }
    }

    // ---- figures per group
    struct FigDef { int g; const char* file; double xmax[3]; };
    const FigDef defs[4] = {{kPi0D, "closure_pi0D", {0.14, 100, 1.25}},
                            {kEtaD, "closure_etaD", {0.55, 180, 1.25}},
                            {kConvPi0, "closure_conv", {0.10, 20, 1.25}},
                            {kModeA, "closure_modeA", {0.60, 180, 1.25}}};
    const char* xt[3] = {"M_{e^{+}e^{-}} [GeV/c^{2}]", "OA_{e^{+}e^{-}} [deg]", "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]"};
    const char* yt[3] = {"d#sigma/dM_{ee} [nb/(GeV/c^{2})]", "d#sigma/dOA [nb/deg]", "d#sigma/dM_{4} [nb/(GeV/c^{2})]"};
    const Audit aud[kNCG] = {{0.894, 0.003, "audit"}, {0.931, 0.005, "audit"}, {1.04, 0.11, "spec"}, {0.952, 0.009, "audit"},
                             {1.06, 0.03, "audit"}, {0.896, 0.004, "spec"}};
    for (const FigDef& d : defs) {
        // the conversion figure shows pi0 (top) and eta (in the ratio row as a second series)
        TCanvas* c = new TCanvas(Form("c_%s", d.file), "", 1800, 1000);
        c->Divide(3, 2, 0.002, 0.002);
        for (int o = 1; o <= 3; ++o) {
            const Bins& bb = *obsB[o];
            const double xmin = (o == 3) ? 0.25 : 0.0, xmax = d.xmax[o - 1];
            std::vector<TH1D*> hs;
            GV A[2], B[2];
            for (int s = 0; s < 2; ++s) {
                A[s] = LoadG(Form("cl_%s_s%d_o%d", kCGKey[d.g], s, o), 0, 0);
                B[s] = LoadG(Form("cl_%s_s%d_o%d", kCGKey[d.g], s, o), 1, 0);
            }
            TVirtualPad* p = c->cd(o);
            PadMargins(p, 0.15, 0.03, 0.10, 0.14);
            p->SetLogy();
            TH1D* hA0 = HistOf(A[0], bb, Form("hA0_%s_%d", d.file, o));
            TH1D* hB0 = HistOf(B[0], bb, Form("hB0_%s_%d", d.file, o));
            TH1D* hA1 = HistOf(A[1], bb, Form("hA1_%s_%d", d.file, o));
            TH1D* hB1 = HistOf(B[1], bb, Form("hB1_%s_%d", d.file, o));
            StyleLine(hA0, kAzure + 2, 3, 1); StyleLine(hB0, kOrange + 7, 2, 1);
            StyleLine(hA1, kAzure + 2, 3, 2); StyleLine(hB1, kOrange + 7, 2, 2);
            hs = {hA0, hB0, hA1, hB1};
            const double ymax = MaxOf(hs, xmin, xmax), ymin = MinPos(hs, xmin, xmax);
            hA0->GetXaxis()->SetRangeUser(xmin, xmax);
            hA0->SetMaximum(ymax * 3); hA0->SetMinimum(std::max(ymin * 0.5, ymax * 1e-5));
            AxisStyle(hA0, xt[o - 1], yt[o - 1]);
            hA0->Draw("hist e");
            hB0->Draw("hist e same");
            hA1->Draw("hist same");
            hB1->Draw("hist same");
            hA0->SetMaximum(ymax * 40);
            if (o == 1) {
                TLegend* lg = new TLegend(0.45, 0.70, 0.97, 0.895);
                lg->SetBorderSize(0); lg->SetFillStyle(1001); lg->SetFillColor(kWhite); lg->SetTextSize(0.040);
                lg->AddEntry(hA0, "std set, pre", "l");
                lg->AddEntry(hB0, "lep set, pre", "l");
                lg->AddEntry(hA1, "std set, pre + OA > 2#circ", "l");
                lg->AddEntry(hB1, "lep set, pre + OA > 2#circ", "l");
                lg->Draw();
            }
            Txt(0.15, 0.925, Form("(%c) %s, e^{+}e^{-}, PDG BR; #sigma = #Sigma w / 221.05 nb^{-1}", 'a' + o - 1,
                                  kCGShort[d.g]), 0.036);
            if (o == 3) Txt(0.18, 0.84, "#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (direct) not simulated", 0.036);
            // ratio
            TVirtualPad* q = c->cd(o + 3);
            PadMargins(q, 0.15, 0.03, 0.06, 0.14);
            TH1D* r0 = RatioHist(B[0], A[0], bb, Form("r0_%s_%d", d.file, o), 10.);
            TH1D* r1 = RatioHist(B[1], A[1], bb, Form("r1_%s_%d", d.file, o), 10.);
            StyleMark(r0, kBlack, 20, 1.2);
            StyleMark(r1, kAzure + 1, 25, 1.2);
            r0->GetXaxis()->SetRangeUser(xmin, xmax);
            r0->SetMinimum(0.5); r0->SetMaximum(1.5);
            AxisStyle(r0, xt[o - 1], "lep / std");
            r0->Draw("e1 p");
            r1->Draw("e1 p same");
            TH1D* r2 = nullptr;
            if (d.g == kConvPi0) {   // eta photon conversions in the same ratio panel
                const GV Ae = LoadG(Form("cl_convEta_s0_o%d", o), 0, 0), Be = LoadG(Form("cl_convEta_s0_o%d", o), 1, 0);
                r2 = RatioHist(Be, Ae, bb, Form("r2_%s_%d", d.file, o), 10.);
                StyleMark(r2, kMagenta + 1, 22, 1.2);
                r2->Draw("e1 p same");
            }
            TLine* l1 = new TLine(xmin, 1, xmax, 1); l1->SetLineStyle(1); l1->SetLineColor(kGray + 2); l1->Draw();
            TLine* l2 = new TLine(xmin, aud[d.g].r, xmax, aud[d.g].r); l2->SetLineStyle(7); l2->SetLineColor(kGreen + 2);
            l2->SetLineWidth(2); l2->Draw();
            if (o == 1) {
                TLegend* lg = new TLegend(0.18, 0.74, 0.97, 0.93);
                lg->SetBorderSize(0); lg->SetFillStyle(1001); lg->SetFillColor(kWhite); lg->SetTextSize(0.038); lg->SetNColumns(2);
                lg->AddEntry(r0, "pre", "pe");
                lg->AddEntry(r1, "pre + OA > 2#circ", "pe");
                if (r2) lg->AddEntry(r2, "#eta-photon conversions, pre", "pe");
                lg->AddEntry(l2, Form("%s total %.3f", aud[d.g].src, aud[d.g].r), "l");
                lg->Draw();
            }
            if (o == 3) Txt(0.18, 0.17, "bins with N_{eff} #geq 10 in both sets; event-group errors", 0.038);
        }
        c->cd(0);
        Save(c, d.file, true);
        delete c;
    }

    // ---- summary figure: total ratios
    {
        struct It { const char* lab; int g; std::vector<std::pair<int, int>> std_, lep_; Audit au; };
        const std::vector<It> items = {
            {"#pi^{0} Dalitz", kPi0D, {{0, 0}}, {{1, 0}}, {0.894, 0.003, ""}},
            {"#eta Dalitz", kEtaD, {{0, 0}}, {{1, 0}}, {0.931, 0.005, ""}},
            {"2 Dalitz e^{+}e^{-}", kDD, {{0, 0}}, {{1, 0}}, {1.19, 0.18, ""}},
            {"2 Dalitz like-sign", kDD, {{0, 1}, {0, 2}}, {{1, 1}, {1, 2}}, {1.09, 0.14, ""}},
            {"conv. #pi^{0} D. #gamma", kConvPi0, {{0, 0}}, {{1, 0}}, {0.952, 0.009, ""}},
            {"conv. #eta D. #gamma", kConvEta, {{0, 0}}, {{1, 0}}, {1.06, 0.03, ""}},
            {"mode A e^{+}e^{-}", kModeA, {{0, 0}}, {{1, 0}}, {0.896, 0.004, ""}},
            {"mode A e^{+}e^{+}", kModeA, {{0, 1}}, {{1, 1}}, {1.10, 0.14, ""}},
            {"mode A e^{-}e^{-}", kModeA, {{0, 2}}, {{1, 2}}, {0.96, 0.08, ""}}};
        const int ni = (int)items.size();
        TCanvas* c = new TCanvas("c_sum", "", 1700, 800);
        PadMargins(c, 0.08, 0.02, 0.10, 0.12);
        TH1D* fr = new TH1D("fr_sum", "", ni, 0, ni);
        fr->SetDirectory(nullptr);
        for (int i = 0; i < ni; ++i) fr->GetXaxis()->SetBinLabel(i + 1, items[i].lab);
        fr->SetMinimum(0.6); fr->SetMaximum(1.45);
        AxisStyle(fr, "", "lep / std (total)");
        fr->GetXaxis()->SetLabelSize(0.040);
        fr->GetXaxis()->LabelsOption("h");
        fr->GetYaxis()->SetTitleOffset(0.7);
        fr->Draw("axis");
        TGraphErrors* gr[3];
        for (int j = 0; j < 3; ++j) gr[j] = new TGraphErrors();
        for (int i = 0; i < ni; ++i) {
            for (int s = 0; s < 2; ++s) {
                const GV a = SumG(Form("cl_%s_s%d_o%d", kCGKey[items[i].g], s, kOTot), items[i].std_);
                const GV b = SumG(Form("cl_%s_s%d_o%d", kCGKey[items[i].g], s, kOTot), items[i].lep_);
                const Rat r = RatioIndep(b.Tot(), b.TotV(), a.Tot(), a.TotV());
                gr[s]->SetPoint(i, i + 0.35 + 0.15 * s, r.r);
                gr[s]->SetPointError(i, 0, r.e);
            }
            gr[2]->SetPoint(i, i + 0.65, items[i].au.r);
            gr[2]->SetPointError(i, 0, items[i].au.e);
        }
        gr[0]->SetMarkerStyle(20); gr[0]->SetMarkerSize(1.5); gr[0]->SetLineWidth(2);
        gr[1]->SetMarkerStyle(25); gr[1]->SetMarkerSize(1.5); gr[1]->SetMarkerColor(kAzure + 1); gr[1]->SetLineColor(kAzure + 1); gr[1]->SetLineWidth(2);
        gr[2]->SetMarkerStyle(33); gr[2]->SetMarkerSize(2.0); gr[2]->SetMarkerColor(kGreen + 2); gr[2]->SetLineColor(kGreen + 2); gr[2]->SetLineWidth(2);
        TLine* l1 = new TLine(0, 1, ni, 1); l1->SetLineColor(kGray + 2); l1->Draw();
        for (int j = 0; j < 3; ++j) gr[j]->Draw("p same");
        TLegend* lg = new TLegend(0.10, 0.72, 0.60, 0.89);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
        lg->AddEntry(gr[0], "FAT production, pre: isBest >= 0, PT3, vz > -160 mm", "pe");
        lg->AddEntry(gr[1], "FAT production, pre + OA > 2 deg", "pe");
        lg->AddEntry(gr[2], "audit (PAT trees, PT3, no vz cut) or spec", "pe");
        lg->Draw();
        Txt(0.08, 0.93, "Lepton set / standard set for the rows the lepton set provides in mode A; PDG BR on both sides; expected 1", 0.036);
        Txt(0.62, 0.85, "pions from the common SMASH final state", 0.032);
        Txt(0.62, 0.80, "(mode A columns: all rows of the class)", 0.032);
        Txt(0.62, 0.75, "event-group errors; simulation only", 0.032);
        Save(c, "closure_summary", true);
        delete c;
    }
}

// ============================================================================
// Part 2: partition mode and w_db
// ============================================================================
void PartPartition(FILE* ft) {
    const char* tcN[5] = {"T1 signal", "T2 physical background", "T3 fake lepton", "T4 different sources", "T1+T2 simulated spectrum"};
    const char* selN[2] = {"pre", "pre+OA>2"};
    fprintf(ft, "# closure_partition: mode B / mode A and w_db / w for e+e- rows, std + lep owned rows added\n");
    fprintf(ft, "# sigma_A = sum w (owned in mode A) / 221.05 nb^-1. Ratios over the same rows with the group-correlated error\n");
    fprintf(ft, "# V(R) = sum_g (b_g - R a_g)^2 / A^2. B/A: weight w; A_db/A: owned A, w_db against w; B_db/A: owned B with w_db against owned A with w\n");
    fprintf(ft, "# A_dbNF, B_dbNF: w_db multiplied by the natural-decay factors NF of src/sim_truth.h (the factors w_db lacks)\n");
    for (int s = 0; s < 2; ++s) {
        fprintf(ft, "\n## selection %s\n", selN[s]);
        fprintf(ft, "%-26s %-14s %20s %10s %18s %18s %18s %18s %18s\n", "template", "M4 window", "sigma_A [nb]", "Neff_A", "B/A",
                "A_db/A", "B_db/A", "A_dbNF/A", "B_dbNF/A");
        for (int tc = 0; tc < 5; ++tc) {
            PV p[5];
            for (int vr = 0; vr < 5; ++vr) p[vr] = LoadP(Form("pa_tc%d_s%d_v%d_b1", tc, s, vr), kStdLepUS);
            for (int i = 0; i <= kNWin; ++i) {
                const bool tot = i == kNWin;
                const double A = p[0].a[i], VA = p[0].aa[i];
                const Rat r0 = RatioPaired(p[0], i), r1 = RatioPaired(p[1], i), r2 = RatioPaired(p[2], i);
                const Rat r3 = RatioPaired(p[3], i), r4 = RatioPaired(p[4], i);
                const std::string wl = tot ? std::string("all M4") : std::string(Form("%.4f-%.4f", kWinEdge[i], kWinEdge[i + 1]));
                fprintf(ft, "%-26s %-14s %20s %10.0f %18s %18s %18s %18s %18s\n", tcN[tc], wl.c_str(),
                        Sig(A / simeng::kLsim, std::sqrt(VA) / simeng::kLsim).c_str(), VA > 0 ? A * A / VA : 0.,
                        PM(r0.r, r0.e).c_str(), PM(r1.r, r1.e).c_str(), PM(r2.r, r2.e).c_str(), PM(r3.r, r3.e).c_str(),
                        PM(r4.r, r4.e).c_str());
            }
        }
        // decomposition of B - A for the physical background
        fprintf(ft, "\n### T2: B - A = (std rows owned only in B) - (lep rows owned only in A), by lepton-pair class, sigma [nb]\n");
        fprintf(ft, "%-14s %14s %14s %14s %14s | %14s %14s %14s %14s | %16s %16s\n", "M4 window", "std+ L1", "std+ L2", "std+ L7",
                "std+ other", "lep- L1", "lep- L2", "lep- L7", "lep- other", "net (B-A)", "std+ / lep- L1");
        const char* lgN[4] = {"L1", "L2", "L7", "other"};
        (void)lgN;
        GV sp[4], lm[4];
        for (int l = 0; l < 4; ++l) {
            sp[l] = LoadG(Form("pd_s%d_l%d_b1", s, l), 0, 0);
            lm[l] = LoadG(Form("pd_s%d_l%d_b1", s, l), 1, 0);
        }
        for (int i = 0; i <= kNWin; ++i) {
            double net = 0;
            fprintf(ft, "%-14s", i == kNWin ? "all M4" : Form("%.4f-%.4f", kWinEdge[i], kWinEdge[i + 1]));
            for (int l = 0; l < 4; ++l) { fprintf(ft, " %14.4g", sp[l].S(i) / simeng::kLsim); net += sp[l].S(i); }
            fprintf(ft, " |");
            for (int l = 0; l < 4; ++l) { fprintf(ft, " %14.4g", lm[l].S(i) / simeng::kLsim); net -= lm[l].S(i); }
            const Rat r = RatioIndep(lm[0].S(i), lm[0].V(i), sp[0].S(i), sp[0].V(i));
            fprintf(ft, " | %16.4g %16s\n", net / simeng::kLsim, r.r > 0 ? PM(1. / r.r, r.e / (r.r * r.r)).c_str() : "-");
        }
    }

    // ---- figures
    for (int s = 0; s < 2; ++s) {
        TCanvas* c = new TCanvas(Form("c_part%d", s), "", 1800, 1100);
        c->Divide(2, 2, 0.002, 0.002);
        const PV pM[3] = {LoadP(Form("pa_tc1_s%d_v0_b0", s), kStdLepUS), LoadP(Form("pa_tc1_s%d_v1_b0", s), kStdLepUS),
                          LoadP(Form("pa_tc1_s%d_v2_b0", s), kStdLepUS)};
        auto mk = [&](const std::vector<double>& x, const std::vector<double>& vx, const char* nm) {
            GV g; g.s = x; g.v = vx; g.n = x;
            return HistOf(g, bM4, nm);
        };
        TH1D* hA = mk(pM[0].a, pM[0].aa, Form("hA_%d", s));
        TH1D* hB = mk(pM[0].b, pM[0].bb, Form("hB_%d", s));
        TH1D* hAd = mk(pM[1].b, pM[1].bb, Form("hAd_%d", s));
        TH1D* hBd = mk(pM[2].b, pM[2].bb, Form("hBd_%d", s));
        StyleLine(hA, kGray + 1, 1, 1, kGray);
        StyleLine(hB, kBlack, 2, 1);
        StyleLine(hAd, kAzure + 2, 2, 2);
        StyleLine(hBd, kOrange + 7, 2, 2);
        const char* xt = "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]";
        for (int pz = 1; pz <= 2; ++pz) {
            TVirtualPad* p = c->cd(pz);
            PadMargins(p, 0.14, 0.03, 0.10, 0.13);
            if (pz == 1) p->SetLogy();
            TH1D* f = (TH1D*)hA->Clone(Form("fA_%d_%d", s, pz));
            f->SetDirectory(nullptr);
            const double xlo = pz == 1 ? 0.30 : 0.40, xhi = pz == 1 ? 1.10 : 0.70;
            const double ym = MaxOf({hA, hB, hAd, hBd}, xlo, xhi);
            f->GetXaxis()->SetRangeUser(xlo, xhi);
            if (pz == 1) { f->SetMaximum(ym * 60); f->SetMinimum(ym * 2e-3); }
            else { f->SetMaximum(ym * 1.15); f->SetMinimum(0); }
            AxisStyle(f, xt, "d#sigma/dM_{4} [nb/(GeV/c^{2})]");
            f->Draw("hist");
            WindowBoxes(pz == 1 ? ym * 2e-3 : 0, pz == 1 ? ym * 60 : ym * 1.15, pz == 1);
            f->Draw("hist same");
            hB->Draw("hist same");
            hAd->Draw("hist same");
            hBd->Draw("hist same");
            gPad->RedrawAxis();
            if (pz == 1) {
                TLegend* lg = new TLegend(0.16, 0.72, 0.97, 0.895);
                lg->SetBorderSize(0); lg->SetFillStyle(1001); lg->SetFillColor(kWhite); lg->SetTextSize(0.038); lg->SetNColumns(2);
                lg->AddEntry(hA, "mode A, w (default)", "f");
                lg->AddEntry(hB, "mode B, w", "l");
                lg->AddEntry(hAd, "mode A, w_{db}", "l");
                lg->AddEntry(hBd, "mode B, w_{db}", "l");
                lg->Draw();
            }
            Txt(0.14, 0.925, Form("(%c) T2 physical background, e^{+}e^{-}, %s; std + lep, #sigma = #Sigma w / 221.05 nb^{-1}",
                                  'a' + pz - 1, s ? "pre + OA > 2#circ" : "pre"), 0.038);
            if (pz == 2) {
                Txt(0.16, 0.84, "windows: bump (orange), peak (green), above peak (blue)", 0.036);
                Txt(0.16, 0.79, "#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (direct) not simulated", 0.036);
            }
        }
        // ratios
        {
            TVirtualPad* p = c->cd(3);
            PadMargins(p, 0.14, 0.03, 0.08, 0.13);
            TH1D* r0 = RatioHistP(pM[0], bM4, Form("rp0_%d", s));
            TH1D* r1 = RatioHistP(pM[1], bM4, Form("rp1_%d", s));
            TH1D* r2 = RatioHistP(pM[2], bM4, Form("rp2_%d", s));
            StyleMark(r0, kBlack, 20, 1.0); StyleMark(r1, kAzure + 2, 24, 1.0); StyleMark(r2, kOrange + 7, 26, 1.0);
            r0->SetMinimum(0.8); r0->SetMaximum(1.75);
            AxisStyle(r0, xt, "variant / (mode A, w)");
            r0->Draw("e1 p");
            WindowBoxes(0.8, 1.75, false);
            r0->Draw("e1 p same"); r1->Draw("e1 p same"); r2->Draw("e1 p same");
            TLine* l1 = new TLine(0.30, 1, 1.10, 1); l1->SetLineColor(kGray + 2); l1->Draw();
            TLegend* lg = new TLegend(0.16, 0.80, 0.97, 0.915);
            lg->SetBorderSize(0); lg->SetFillStyle(1001); lg->SetFillColor(kWhite); lg->SetTextSize(0.040); lg->SetNColumns(3);
            lg->AddEntry(r0, "B, w", "pe"); lg->AddEntry(r1, "A, w_{db}", "pe"); lg->AddEntry(r2, "B, w_{db}", "pe");
            lg->Draw();
            Txt(0.14, 0.935, "(c) ratio to the default per 10 MeV bin, group-correlated errors", 0.038);
            gPad->RedrawAxis();
        }
        // decomposition of B - A
        {
            TVirtualPad* p = c->cd(4);
            PadMargins(p, 0.14, 0.03, 0.08, 0.13);
            const int col[4] = {kOrange + 7, kMagenta + 1, kGreen + 2, kAzure + 2};
            const char* nm[4] = {"L1 #pi^{0} Dalitz", "L2 #eta Dalitz", "L7 two Dalitz", "other"};
            TH1D* hp[4]; TH1D* hm[4];
            TH1D* net = nullptr;
            double ym = 0;
            for (int l = 0; l < 4; ++l) {
                GV a = LoadG(Form("pd_s%d_l%d_b0", s, l), 0, 0), b = LoadG(Form("pd_s%d_l%d_b0", s, l), 1, 0);
                hp[l] = HistOf(a, bM4, Form("hp_%d_%d", s, l));
                hm[l] = HistOf(b, bM4, Form("hm_%d_%d", s, l));
                hm[l]->Scale(-1);
                StyleLine(hp[l], col[l], 2, 1); StyleLine(hm[l], col[l], 2, 2);
                if (!net) { net = (TH1D*)hp[l]->Clone(Form("net_%d", s)); net->SetDirectory(nullptr); net->Add(hm[l]); }
                else { net->Add(hp[l]); net->Add(hm[l]); }
                ym = std::max(ym, std::max(hp[l]->GetMaximum(), -hm[l]->GetMinimum()));
            }
            StyleMark(net, kBlack, 20, 0.9);
            hp[0]->SetMaximum(ym * 1.25); hp[0]->SetMinimum(-ym * 1.25);
            AxisStyle(hp[0], xt, "d#sigma/dM_{4} [nb/(GeV/c^{2})]");
            hp[0]->Draw("hist");
            WindowBoxes(-ym * 1.25, ym * 1.25, false);
            for (int l = 0; l < 4; ++l) { hp[l]->Draw("hist same"); hm[l]->Draw("hist same"); }
            net->Draw("e1 p same");
            TLine* l0 = new TLine(0.30, 0, 1.10, 0); l0->SetLineColor(kGray + 2); l0->Draw();
            hp[0]->SetMaximum(ym * 1.9);
            TLegend* lg = new TLegend(0.16, 0.73, 0.97, 0.915);
            lg->SetBorderSize(0); lg->SetFillStyle(1001); lg->SetFillColor(kWhite); lg->SetTextSize(0.036); lg->SetNColumns(3);
            for (int l = 0; l < 4; ++l) lg->AddEntry(hp[l], nm[l], "l");
            lg->AddEntry(net, "net B - A", "pe");
            lg->Draw();
            Txt(0.14, 0.935, "(d) B - A of T2: + std rows owned only in B (solid), - lep rows owned only in A (dashed)", 0.036);
            gPad->RedrawAxis();
        }
        c->cd(0);
        Save(c, s ? "closure_partition_oa2" : "closure_partition_pre", true);
        delete c;
    }
}

// ============================================================================
// Part 3: flags and rules-v1 gaps
// ============================================================================
void PartFlags(FILE* ft) {
    const char* tpN[4] = {"T1", "T2", "T3", "T4"};
    const char* ocN[3] = {"owned A", "owned B only", "owned by none"};
    const char* hypN[3] = {"e+e-", "e+e+", "e-e-"};
    const char* bitN[15] = {"0 weight mismatch", "1 unknown forced (reserved)", "2 chain factor not applied",
                            "3 unknown code", "4 ambiguous code", "5 e+- in a pion slot", "6 second-decay product",
                            "7 several forced decays", "8 weight without forced code", "9 code of another set",
                            "10 ghost lepton", "11 ghost pion", "12 same particle in both lepton slots",
                            "13 rho0/omega/phi + FC (no set)", "14 v1 code gap (phi K0K0bar / phi->eta gamma)"};
    const int winShow[4] = {1, 3, 5, kNWin};   // bump, peak, above peak, all
    fprintf(ft, "# closure_flags: rows and weight per sim_flags bit after pre (and pre+OA>2), per template and ownership\n");
    fprintf(ft, "# sigma = sum w / 221.05 nb^-1 (w = sim_w as stored, PDG BR). frac = sigma / sigma of the same template owned in mode A\n");
    fprintf(ft, "# windows: bump 0.455-0.5225, peak 0.5336-0.5542, above peak 0.58-0.70 GeV/c2; all = every M4\n");
    for (int s = 0; s < 2; ++s) {
        fprintf(ft, "\n## selection %s\n", s ? "pre+OA>2" : "pre");
        for (int h = 0; h < 3; ++h) {
            fprintf(ft, "\n### %s (std + lep)\n", hypN[h]);
            fprintf(ft, "%-44s %-4s %-13s", "bit", "tmpl", "ownership");
            for (int wi : winShow) fprintf(ft, " | %-34s", wi == kNWin ? "all M4: rows, sigma [nb], frac" : Form("%s: rows, sigma, frac", kWinName[wi]));
            fprintf(ft, "\n");
            GV den[4];
            for (int tp = 0; tp < 4; ++tp) den[tp] = SumG(Form("tt_s%d_tp%d_o0", s, tp), {{0, h}, {1, h}});
            for (int bit = 0; bit < 15; ++bit)
                for (int tp = 0; tp < 4; ++tp)
                    for (int oc = 0; oc < 3; ++oc) {
                        const GV g = SumG(Form("fl_s%d_b%d_tp%d_o%d", s, bit, tp, oc), {{0, h}, {1, h}});
                        if (g.empty() || g.TotN() == 0) continue;
                        fprintf(ft, "%-44s %-4s %-13s", bitN[bit], tpN[tp], ocN[oc]);
                        for (int wi : winShow) {
                            const double fr = den[tp].S(wi) > 0 ? g.S(wi) / den[tp].S(wi) : 0;
                            fprintf(ft, " | %7.0f %14s %9.2e", g.N(wi), Sig(g.S(wi) / simeng::kLsim, std::sqrt(g.V(wi)) / simeng::kLsim).c_str(), fr);
                        }
                        fprintf(ft, "\n");
                    }
        }
    }
    // ---- gap effects on the templates (e+e-)
    const char* efN[6] = {"bit 13: rows owned by no set (missing weight)", "bit 14: std rows dropped in A, no set provides",
                          "bit 14: phi->eta gamma, eta factor missing (W'-W)", "bit 2: chain factor not applied (upper bound = all)",
                          "w_db without natural factors: w_db (NF-1)", "bit 12: same particle twice (owned A)"};
    fprintf(ft, "\n\n# gap effects on the e+e- templates: change / template sum (owned in mode A, w; for the w_db gap: / sum w_db)\n");
    fprintf(ft, "# a positive number = weight that is missing (or would be added) if the gap were closed\n");
    for (int s = 0; s < 2; ++s) {
        fprintf(ft, "\n## selection %s\n", s ? "pre+OA>2" : "pre");
        fprintf(ft, "%-52s %-4s", "gap", "tmpl");
        for (int i = 0; i <= kNWin; ++i) fprintf(ft, " %22s", i == kNWin ? "all M4" : kWinName[i]);
        fprintf(ft, "\n");
        for (int e = 0; e < 6; ++e)
            for (int tp = 0; tp < 4; ++tp) {
                const GV g = SumG(Form("ef_s%d_e%d_tp%d_b1", s, e, tp), kStdLepUS);
                const GV d = e == 4 ? SumG(Form("wdb_s%d_tp%d_b1", s, tp), kStdLepUS) : SumG(Form("tt_s%d_tp%d_o0", s, tp), kStdLepUS);
                fprintf(ft, "%-52s %-4s", efN[e], tpN[tp]);
                for (int i = 0; i <= kNWin; ++i) {
                    const double f = d.S(i) > 0 ? g.S(i) / d.S(i) : 0, fe = d.S(i) > 0 ? std::sqrt(g.V(i)) / d.S(i) : 0;
                    fprintf(ft, " %22s", PM(100 * f, 100 * fe, "%.3g%%").c_str());
                }
                fprintf(ft, "\n");
            }
    }
    // like-sign: bit 12 share of the owned weight
    fprintf(ft, "\n# like-sign trees: bit 12 rows (clone pairs) among the owned-A rows after pre, all M4\n");
    for (int h = 1; h < 3; ++h) {
        double allW = 0, cl = 0, allN = 0, clN = 0;
        for (int tp = 0; tp < 4; ++tp) {
            const GV d = SumG(Form("tt_s0_tp%d_o0", tp), {{0, h}, {1, h}});
            const GV g = SumG(Form("ef_s0_e5_tp%d_b1", tp), {{0, h}, {1, h}});
            allW += d.Tot(); allN += d.TotN(); cl += g.Tot(); clN += g.TotN();
        }
        fprintf(ft, "%s: %.0f of %.0f rows, sigma %.4g of %.4g nb (%.2f %%)\n", hypN[h], clN, allN, cl / simeng::kLsim,
                allW / simeng::kLsim, allW > 0 ? 100 * cl / allW : 0);
    }

    // ---- figure: effects per M4 bin relative to T1+T2
    TCanvas* c = new TCanvas("c_flags", "", 1500, 1100);
    TPad* p1 = new TPad("p1", "", 0, 0.42, 1, 1);
    TPad* p2 = new TPad("p2", "", 0, 0, 1, 0.42);
    p1->Draw(); p2->Draw();
    PadMargins(p1, 0.13, 0.03, 0.08, 0.02);
    PadMargins(p2, 0.13, 0.03, 0.03, 0.22);
    GV sim;   // T1 + T2 owned A, w
    {
        const PV p = LoadP("pa_tc4_s0_v0_b0", kStdLepUS);
        sim.s = p.a; sim.v = p.aa; sim.n = p.n;
    }
    GV wdbT12 = SumG("wdb_s0_tp0_b0", kStdLepUS);
    {
        GV x = SumG("wdb_s0_tp1_b0", kStdLepUS);
        AddVec(wdbT12.s, x.s); AddVec(wdbT12.v, x.v); AddVec(wdbT12.n, x.n);
    }
    TH1D* hs = HistOf(sim, bM4, "h_sim_fl");
    StyleLine(hs, kGray + 1, 1, 1, kGray);
    const int ecol[6] = {kMagenta + 1, kGreen + 2, kGreen + 3, kOrange + 7, kAzure + 2, kRed + 1};
    const char* elab[6] = {"bit 13: #rho^{0}/#omega/#phi + FC, owned by no set", "bit 14: std rows dropped, no set provides",
                           "bit 14: #phi#rightarrow#eta#gamma, |W'-W|", "bit 2: chain factor (all such rows)",
                           "w_{db} natural factors: |w_{db}(NF-1)|", "bit 12: same particle twice"};
    TH1D* he[6];
    TH1D* hr[6];
    for (int e = 0; e < 6; ++e) {
        GV g;
        for (int tp = 0; tp < 2; ++tp) {   // effects on T1 + T2
            const GV x = SumG(Form("ef_s0_e%d_tp%d_b0", e, tp), kStdLepUS);
            AddVec(g.s, x.s); AddVec(g.v, x.v); AddVec(g.n, x.n);
        }
        if (g.empty()) { g.s.assign(bM4.n() + 1, 0.); g.v = g.s; g.n = g.s; }
        for (auto& x : g.s) x = std::fabs(x);
        he[e] = HistOf(g, bM4, Form("he_%d", e));
        StyleLine(he[e], ecol[e], 2, e == 4 ? 2 : 1);
        hr[e] = new TH1D(Form("hr_%d", e), "", bM4.n(), bM4.e.data());
        hr[e]->SetDirectory(nullptr);
        const GV& dd = e == 4 ? wdbT12 : sim;
        for (int i = 0; i < bM4.n(); ++i) {
            const double v = dd.S(i) > 0 ? 100. * g.S(i) / dd.S(i) : 0;
            hr[e]->SetBinContent(i + 1, v > 0 ? v : 1e-9);
            hr[e]->SetBinError(i + 1, dd.S(i) > 0 ? 100. * std::sqrt(g.V(i)) / dd.S(i) : 0);
        }
        StyleMark(hr[e], ecol[e], 20 + e, 0.9);
    }
    p1->cd(); p1->SetLogy();
    const double ym = hs->GetMaximum();
    hs->SetMaximum(ym * 3000); hs->SetMinimum(ym * 1e-7);
    AxisStyle(hs, "", "d#sigma/dM_{4} [nb/(GeV/c^{2})]", 0.75);
    hs->GetYaxis()->SetTitleOffset(1.1);
    hs->GetXaxis()->SetLabelSize(0);
    hs->Draw("hist");
    WindowBoxes(ym * 1e-7, ym * 3000, true);
    hs->Draw("hist same");
    for (int e = 0; e < 6; ++e) he[e]->Draw("hist same");
    gPad->RedrawAxis();
    TLegend* lg = new TLegend(0.14, 0.70, 0.97, 0.905);
    lg->SetBorderSize(0); lg->SetFillStyle(1001); lg->SetFillColor(kWhite); lg->SetTextSize(0.032); lg->SetNColumns(2);
    lg->AddEntry(hs, "T1 + T2 (simulated spectrum), mode A, w", "f");
    for (int e = 0; e < 6; ++e) lg->AddEntry(he[e], elab[e], "l");
    lg->Draw();
    Txt(0.13, 0.935, "Open rules-v1 gaps, e^{+}e^{-}, pre; std + lep owned rows (mode A), #sigma = #Sigma w / 221.05 nb^{-1}; #eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (direct) not simulated", 0.032);
    p2->cd(); p2->SetLogy();
    TH1D* fr = (TH1D*)hr[0]->Clone("fr_flags");
    fr->SetDirectory(nullptr);
    fr->Reset();
    fr->SetMinimum(1e-4); fr->SetMaximum(300);
    AxisStyle(fr, "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]", "|change| / template [%]", 1.3);
    fr->GetYaxis()->SetTitleOffset(0.5);
    fr->Draw("axis");
    WindowBoxes(1e-4, 300, true);
    for (int e = 0; e < 6; ++e) hr[e]->Draw("p same");
    TLine* l1 = new TLine(0.30, 1, 1.10, 1); l1->SetLineStyle(2); l1->SetLineColor(kBlack); l1->SetLineWidth(2); l1->Draw();
    Txt(0.12, 0.88, "relative to T1+T2 (w_{db} gap: relative to T1+T2 with w_{db}); dashed line: 1 %", 0.055);
    gPad->RedrawAxis();
    c->cd();
    Save(c, "closure_flags", true);
    delete c;
}

// ============================================================================
// Part 4: ghost leptons in the signal classes
// ============================================================================
struct SigRow { int cls, gl, oa2, grp, set; float m4, m4t; double w; };
struct Pk { double mean = 0, emean = 0, sig = 0, esig = 0, med = 0, hw68 = 0, chi2 = 0; int ndf = 0; };
Pk PeakOf(const std::vector<SigRow>& v, const char* name0) {
    Pk r;
    const std::string name(name0);
    TH1D h(name.c_str(), "", 45, 0.50, 0.59);
    h.Sumw2();
    std::vector<std::pair<double, double>> mw;
    for (const auto& x : v) { h.Fill(x.m4, x.w); if (x.m4 > 0.50 && x.m4 < 0.59) mw.push_back({x.m4, x.w}); }
    double mu = 0.5455, sg = 0.006;
    TF1 f(("f_" + name).c_str(), "gaus", 0.50, 0.59);
    for (int it = 0; it < 5; ++it) {
        f.SetParameters(h.GetMaximum(), mu, sg);
        h.Fit(&f, "QS0", "", mu - 2 * sg, mu + 2 * sg);
        mu = f.GetParameter(1); sg = std::fabs(f.GetParameter(2));
    }
    r.mean = mu; r.emean = f.GetParError(1); r.sig = sg; r.esig = f.GetParError(2);
    r.chi2 = f.GetChisquare(); r.ndf = f.GetNDF();
    std::sort(mw.begin(), mw.end());
    double tot = 0;
    for (auto& x : mw) tot += x.second;
    auto q = [&](double p) { double c = 0; for (auto& x : mw) { c += x.second; if (c >= p * tot) return x.first; } return mw.empty() ? 0. : mw.back().first; };
    r.med = q(0.5); r.hw68 = 0.5 * (q(0.84) - q(0.16));
    return r;
}
struct WS { double S = 0, V = 0, N = 0; };
WS WinSum(const std::vector<SigRow>& v, double lo, double hi) {
    WS r;
    std::map<std::pair<int, int>, double> g;
    for (const auto& x : v) if (x.m4 >= lo && x.m4 < hi) { r.S += x.w; r.N += 1; g[{x.set, x.grp}] += x.w; }
    for (auto& kv : g) r.V += kv.second * kv.second;
    return r;
}
void PartGhost(FILE* ft) {
    std::vector<SigRow> rows;
    for (int k = 0; k < 6; ++k) {
        TTree* t = gF[k]->Get<TTree>("sig");
        if (!t) continue;
        Int_t cls, gl, oa2, grp, set; Float_t m4, m4t; Double_t w;
        t->SetBranchAddress("cls", &cls); t->SetBranchAddress("gl", &gl); t->SetBranchAddress("oa2", &oa2);
        t->SetBranchAddress("grp", &grp); t->SetBranchAddress("set", &set); t->SetBranchAddress("m4", &m4);
        t->SetBranchAddress("m4t", &m4t); t->SetBranchAddress("w", &w);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); rows.push_back({cls, gl, oa2, grp + 10000000 * k, set, m4, m4t, w}); }
    }
    const char* clsN[5] = {"", "", "C2 bump, pi0 Dalitz", "C3 bump, conversion", "C4 peak, conversion"};
    fprintf(ft, "# closure_ghost: signal rows (template 1) with a ghost track in a lepton slot (sim_flags bit 10) against clean rows\n");
    fprintf(ft, "# e+e-, std set (the lepton set has no signal rows), sigma = sum w / 221.05 nb^-1, event-group errors\n");
    fprintf(ft, "# windows: bump 0.455-0.5225, peak 0.5336-0.5542 GeV/c2\n");
    for (int s = 0; s < 2; ++s) {
        fprintf(ft, "\n## selection %s\n", s ? "pre+OA>2" : "pre");
        fprintf(ft, "%-22s %-6s %7s %18s %18s %18s %8s\n", "class", "rows", "n", "sigma all [nb]", "sigma bump win", "sigma peak win", "ghost%");
        for (int c = 2; c <= 4; ++c)
            for (int gsel = 0; gsel < 3; ++gsel) {   // 0 all, 1 clean, 2 ghost
                std::vector<SigRow> v;
                for (auto& x : rows) if (x.cls == c && (s == 0 || x.oa2) && (gsel == 0 || (gsel == 1 ? !x.gl : x.gl))) v.push_back(x);
                const WS a = WinSum(v, -1e9, 1e9), b = WinSum(v, kWinEdge[1], kWinEdge[2]), p = WinSum(v, kWinEdge[3], kWinEdge[4]);
                std::vector<SigRow> va;
                for (auto& x : rows) if (x.cls == c && (s == 0 || x.oa2)) va.push_back(x);
                const WS all = WinSum(va, -1e9, 1e9);
                fprintf(ft, "%-22s %-6s %7.0f %18s %18s %18s %8.1f\n", clsN[c], gsel == 0 ? "all" : gsel == 1 ? "clean" : "ghost", a.N,
                        Sig(a.S / simeng::kLsim, std::sqrt(a.V) / simeng::kLsim).c_str(),
                        Sig(b.S / simeng::kLsim, std::sqrt(b.V) / simeng::kLsim).c_str(),
                        Sig(p.S / simeng::kLsim, std::sqrt(p.V) / simeng::kLsim).c_str(), all.S > 0 ? 100 * a.S / all.S : 0);
            }
        // peak of C4 and bump shape of C2 + C3
        fprintf(ft, "\n%-34s %22s %22s %12s %12s %10s\n", "C4 peak, conversion", "Gauss mean [GeV/c2]", "Gauss sigma [MeV/c2]",
                "median", "hw68 [MeV]", "chi2/ndf");
        for (int gsel = 0; gsel < 3; ++gsel) {
            std::vector<SigRow> v;
            for (auto& x : rows) if (x.cls == 4 && (s == 0 || x.oa2) && (gsel == 0 || (gsel == 1 ? !x.gl : x.gl))) v.push_back(x);
            const Pk p = PeakOf(v, Form("pk_%d_%d", s, gsel));
            fprintf(ft, "%-34s %22s %22s %12.5f %12.2f %6.1f/%d\n", gsel == 0 ? "all" : gsel == 1 ? "clean" : "ghost",
                    PM(p.mean, p.emean, "%.5f").c_str(), PM(1000 * p.sig, 1000 * p.esig, "%.2f").c_str(), p.med, 1000 * p.hw68, p.chi2, p.ndf);
        }
        fprintf(ft, "\n%-34s %14s %14s %14s %14s\n", "M4 - M4_true (all signal classes)", "mean [MeV]", "rms [MeV]",
                "|d|>10 MeV %", "rows");
        for (int gsel = 1; gsel < 3; ++gsel) {
            double sw = 0, sx = 0, sxx = 0, out = 0, n = 0;
            for (auto& x : rows) {
                if ((s == 1 && !x.oa2) || (gsel == 1 ? x.gl : !x.gl) || x.m4t < 0) continue;
                const double d = 1000. * (x.m4 - x.m4t);
                n += 1;
                if (std::fabs(d) > 10) { out += x.w; }
                sw += x.w;
                if (std::fabs(d) < 40) { sx += x.w * d; sxx += x.w * d * d; }
            }
            double sw40 = 0;
            for (auto& x : rows) if (!((s == 1 && !x.oa2) || (gsel == 1 ? x.gl : !x.gl) || x.m4t < 0) && std::fabs(1000. * (x.m4 - x.m4t)) < 40) sw40 += x.w;
            const double mean = sw40 > 0 ? sx / sw40 : 0, rms = sw40 > 0 ? std::sqrt(std::max(0., sxx / sw40 - mean * mean)) : 0;
            fprintf(ft, "%-34s %14.2f %14.2f %14.1f %14.0f\n", gsel == 1 ? "clean" : "ghost", mean, rms, sw > 0 ? 100 * out / sw : 0, n);
        }
        fprintf(ft, "# (mean and rms inside |d| < 40 MeV)\n");
    }

    // ---- figure
    TCanvas* c = new TCanvas("c_ghost", "", 1800, 650);
    c->Divide(3, 1, 0.002, 0.002);
    auto H = [&](const char* nm, int cl1, int cl2, int gsel, int nb, double lo, double hi, bool res) {
        TH1D* h = new TH1D(nm, "", nb, lo, hi);
        h->SetDirectory(nullptr); h->Sumw2();
        for (auto& x : rows) {
            if (x.cls != cl1 && x.cls != cl2) continue;
            if (gsel == 1 && x.gl) continue;
            if (gsel == 2 && !x.gl) continue;
            if (res) { if (x.m4t > 0) h->Fill(1000. * (x.m4 - x.m4t), x.w); }
            else h->Fill(x.m4, x.w);
        }
        if (!res) h->Scale(1. / (simeng::kLsim * h->GetBinWidth(1)));
        return h;
    };
    {
        TVirtualPad* p = c->cd(1);
        PadMargins(p, 0.16, 0.03, 0.10, 0.14);
        TH1D* a = H("g_c4a", 4, 4, 0, 45, 0.50, 0.59, false);
        TH1D* b = H("g_c4c", 4, 4, 1, 45, 0.50, 0.59, false);
        TH1D* g = H("g_c4g", 4, 4, 2, 45, 0.50, 0.59, false);
        StyleLine(a, kBlack, 2, 1); StyleLine(b, kGreen + 2, 1, 1, kGreen - 9); StyleLine(g, kMagenta + 1, 3, 2);
        a->SetMaximum(a->GetMaximum() * 1.5); a->SetMinimum(0);
        AxisStyle(a, "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]", "d#sigma/dM_{4} [nb/(GeV/c^{2})]");
        a->GetXaxis()->SetNdivisions(505);
        a->Draw("hist e"); b->Draw("hist same"); a->Draw("hist e same"); g->Draw("hist e same");
        TLegend* lg = new TLegend(0.18, 0.62, 0.62, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.045);
        lg->AddEntry(a, "all rows", "l"); lg->AddEntry(b, "clean leptons", "f"); lg->AddEntry(g, "ghost lepton", "l");
        lg->Draw();
        Txt(0.16, 0.925, "(a) C4 peak (#eta#rightarrow#pi^{+}#pi^{-}#gamma, #gamma#rightarrowe^{+}e^{-}), pre, std set", 0.040);
        Txt(0.64, 0.84, "#sigma = #Sigma w / 221.05 nb^{-1}", 0.036);
        Txt(0.64, 0.79, "#eta#rightarrow#pi^{+}#pi^{-}e^{+}e^{-} (direct)", 0.036);
        Txt(0.64, 0.75, "not simulated", 0.036);
    }
    {
        TVirtualPad* p = c->cd(2);
        PadMargins(p, 0.16, 0.03, 0.10, 0.14);
        TH1D* a2 = H("g_c2c", 2, 2, 1, 36, 0.40, 0.58, false);
        TH1D* g2 = H("g_c2g", 2, 2, 2, 36, 0.40, 0.58, false);
        TH1D* a3 = H("g_c3c", 3, 3, 1, 36, 0.40, 0.58, false);
        TH1D* g3 = H("g_c3g", 3, 3, 2, 36, 0.40, 0.58, false);
        StyleLine(a2, kOrange + 7, 3, 1); StyleLine(g2, kOrange + 7, 2, 2);
        StyleLine(a3, kMagenta + 1, 3, 1); StyleLine(g3, kMagenta + 1, 2, 2);
        a2->SetMaximum(std::max(a2->GetMaximum(), a3->GetMaximum()) * 1.45); a2->SetMinimum(0);
        AxisStyle(a2, "M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}]", "d#sigma/dM_{4} [nb/(GeV/c^{2})]");
        a2->GetXaxis()->SetNdivisions(505);
        a2->Draw("hist e"); g2->Draw("hist e same"); a3->Draw("hist e same"); g3->Draw("hist e same");
        TLegend* lg = new TLegend(0.18, 0.62, 0.70, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
        lg->AddEntry(a2, "C2 bump #pi^{0} Dalitz, clean", "l"); lg->AddEntry(g2, "C2, ghost lepton", "l");
        lg->AddEntry(a3, "C3 bump conversion, clean", "l"); lg->AddEntry(g3, "C3, ghost lepton", "l");
        lg->Draw();
        Txt(0.16, 0.925, "(b) bump classes (#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}), pre, std set", 0.040);
    }
    {
        TVirtualPad* p = c->cd(3);
        PadMargins(p, 0.16, 0.03, 0.10, 0.14);
        p->SetLogy();
        TH1D* rc = H("g_rc", 2, 4, 1, 40, -40, 40, true);   // classes 2..4 (filter below)
        TH1D* rg = H("g_rg", 2, 4, 2, 40, -40, 40, true);
        rc->Reset(); rg->Reset();
        for (auto& x : rows) {
            if (x.cls < 2 || x.cls > 4 || x.m4t <= 0) continue;
            (x.gl ? rg : rc)->Fill(1000. * (x.m4 - x.m4t), x.w);
        }
        for (TH1D* h : {rc, rg}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
        StyleLine(rc, kGreen + 2, 3, 1, 0); StyleLine(rg, kMagenta + 1, 3, 2);
        rc->SetMaximum(1.0); rc->SetMinimum(3e-4);
        AxisStyle(rc, "M_{4} - M_{4}^{true} [MeV/c^{2}]", "fraction per 2 MeV/c^{2}");
        rc->Draw("hist e"); rg->Draw("hist e same");
        TLegend* lg = new TLegend(0.18, 0.74, 0.97, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
        lg->AddEntry(rc, "signal rows C2-C4, clean leptons", "l"); lg->AddEntry(rg, "signal rows C2-C4, ghost lepton", "l");
        lg->Draw();
        Txt(0.16, 0.925, "(c) M4 resolution, pre, shapes normalised to 1", 0.040);
    }
    c->cd(0);
    Save(c, "closure_ghost", true);
    delete c;
}

// ============================================================================
// Part 5: float-code dictionary
// ============================================================================
void PartCodes(FILE* ft) {
    struct Cnt { long long n[6][2] = {{0}}; long long tot = 0; };
    std::map<float, Cnt> vals[3];
    for (int k = 0; k < 6; ++k) {
        TTree* t = gF[k]->Get<TTree>("codes");
        Int_t col, styp; Float_t val; Long64_t cnt;
        t->SetBranchAddress("col", &col); t->SetBranchAddress("styp", &styp); t->SetBranchAddress("val", &val);
        t->SetBranchAddress("cnt", &cnt);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            Cnt& c = vals[col][val];
            c.n[k][styp] += cnt; c.tot += cnt;
        }
    }
    const char* colN[3] = {"geninfo", "geninfo1", "geninfo2"};
    const char* trN[6] = {"std_us", "std_pp", "std_mm", "lep_us", "lep_pp", "lep_mm"};
    fprintf(ft, "# closure_codes: every stored generator-code value of the four slots in the six FAT outputs (all rows, no selection)\n");
    fprintf(ft, "# mapped with simtruth::CodeCandidates (src/sim_truth.h): candidates = exact integer codes of the dictionary that round to the stored float\n");
    fprintf(ft, "# counts = track slots (4 per row); per tree: lepton slots / pion slots\n\n");
    long long nUnk[3] = {0}, nAmb[3] = {0}, vUnk[3] = {0}, vAmb[3] = {0};
    std::map<int, std::pair<long long, long long>> famStat[3];   // family -> (distinct values, slots)
    std::map<int, std::pair<double, double>> famRange[3];
    for (int col = 0; col < 3; ++col)
        for (auto& kv : vals[col]) {
            simtruth::Code c[8];
            const int n = simtruth::CodeCandidates(col, kv.first, c, 8);
            if (n == 0) { vUnk[col]++; nUnk[col] += kv.second.tot; }
            if (n > 1) { vAmb[col]++; nAmb[col] += kv.second.tot; }
            const int fam = n >= 1 ? c[0].family : 0;
            famStat[col][fam].first++; famStat[col][fam].second += kv.second.tot;
            auto& rg = famRange[col][fam];
            if (famStat[col][fam].first == 1) rg = {kv.first, kv.first};
            rg.first = std::min(rg.first, (double)kv.first); rg.second = std::max(rg.second, (double)kv.first);
        }
    fprintf(ft, "## summary per column\n%-9s %10s %14s %12s %14s %12s %14s\n", "column", "distinct", "slots", "unknown val", "unknown slots",
            "ambig. val", "ambig. slots");
    for (int col = 0; col < 3; ++col) {
        long long tot = 0;
        for (auto& kv : vals[col]) tot += kv.second.tot;
        fprintf(ft, "%-9s %10zu %14lld %12lld %14lld %12lld %14lld\n", colN[col], vals[col].size(), tot, vUnk[col], nUnk[col], vAmb[col], nAmb[col]);
    }
    fprintf(ft, "\n## families per column: distinct values, slots, smallest and largest stored value\n");
    for (int col = 0; col < 3; ++col)
        for (auto& kv : famStat[col])
            fprintf(ft, "%-9s %-36s %10lld %14lld %14.0f %14.0f\n", colN[col], simtruth::CodeFamilyName(kv.first), kv.second.first,
                    kv.second.second, famRange[col][kv.first].first, famRange[col][kv.first].second);
    // full list for geninfo and geninfo1; geninfo2: only values that are not plain counters / process numbers
    for (int col = 0; col < 3; ++col) {
        fprintf(ft, "\n## %s: %s\n", colN[col], col < 2 ? "every distinct value" : "values outside the families 'counter' and '1000000 + SMASH process' (none expected)");
        fprintf(ft, "%14s %16s %-34s %-52s", "stored", "candidate(s)", "family", "dictionary text");
        for (int k = 0; k < 6; ++k) fprintf(ft, " %17s", trN[k]);
        fprintf(ft, "\n");
        for (auto& kv : vals[col]) {
            simtruth::Code c[8];
            const int n = simtruth::CodeCandidates(col, kv.first, c, 8);
            const int fam = n >= 1 ? c[0].family : 0;
            if (col == 2 && n == 1 && (fam == simtruth::kFamCounter || fam == simtruth::kFamSmashProcess)) continue;
            std::string cs;
            for (int j = 0; j < std::min(n, 8); ++j) cs += Form("%s%lld", j ? "," : "", c[j].code);
            if (n == 0) cs = "UNKNOWN";
            const char* txt = "";
            for (const auto& d : simtruth::kDictCodes) if (n >= 1 && d.column == col && d.code == c[0].code) txt = d.text;
            std::string fname = n == 0 ? "-" : simtruth::CodeFamilyName(fam);
            if (n > 1) fname += " (AMBIGUOUS)";
            if (col == 0 && fam == simtruth::kFamSmashParticle && n >= 1)
                txt = (c[0].code == 1016 || c[0].code == 1610) ? "SMASH K0 / K0bar (code 1016 / 1610)"
                                                                 : Form("SMASH particle, GEANT id %lld", c[0].code / 101);
            fprintf(ft, "%14.0f %16s %-34s %-52s", kv.first, cs.c_str(), fname.c_str(), txt);
            for (int k = 0; k < 6; ++k) fprintf(ft, " %8lld/%-8lld", kv.second.n[k][0], kv.second.n[k][1]);
            fprintf(ft, "\n");
        }
    }
    // cross-check with sim_flags bits 3 and 4
    fprintf(ft, "\n## cross-check per tree with sim_flags (all rows)\n%-8s %10s %14s %14s %14s %14s %14s %14s\n", "tree", "rows",
            "rows unknown", "rows ambig.", "bit 3 rows", "bit 4 rows", "bit3 mismatch", "bit4 mismatch");
    for (int k = 0; k < 6; ++k) {
        TObjString* os = gF[k]->Get<TObjString>("counters");
        std::map<std::string, double> m;
        std::string sx = os->GetString().Data();
        size_t p0 = 0;
        while (p0 < sx.size()) {
            size_t e = sx.find('\n', p0); if (e == std::string::npos) e = sx.size();
            const std::string line = sx.substr(p0, e - p0);
            const size_t tb = line.find('\t');
            if (tb != std::string::npos) m[line.substr(0, tb)] = atof(line.substr(tb + 1).c_str());
            p0 = e + 1;
        }
        fprintf(ft, "%-8s %10.0f %14.0f %14.0f %14.0f %14.0f %14.0f %14.0f\n", trN[k], m["rows"], m["code_rows_unknown"],
                m["code_rows_ambiguous"], m["allrows_bit3"], m["allrows_bit4"], m["code_bit3_mismatch"], m["code_bit4_mismatch"]);
    }
}

// counters of all pass files in one table
void PartCounters(FILE* ft) {
    std::map<std::string, std::vector<double>> all;
    for (int k = 0; k < 6; ++k) {
        TObjString* os = gF[k]->Get<TObjString>("counters");
        std::string sx = os->GetString().Data();
        size_t p0 = 0;
        while (p0 < sx.size()) {
            size_t e = sx.find('\n', p0); if (e == std::string::npos) e = sx.size();
            const std::string line = sx.substr(p0, e - p0);
            const size_t tb = line.find('\t');
            if (tb != std::string::npos) {
                auto& v = all[line.substr(0, tb)];
                if (v.empty()) v.assign(6, 0.);
                v[k] = atof(line.substr(tb + 1).c_str());
            }
            p0 = e + 1;
        }
    }
    fprintf(ft, "# closure_counters: row counters of the pass (all rows / after pre); cl_<group>_own<A><B>: ownership of the closure rows;\n");
    fprintf(ft, "# bit2_*: rows with bit 2 after pre (own0/own1 = owned in A no/yes), samekey = omega lepton and second-decay product share (geninfo, geninfo2)\n");
    fprintf(ft, "%-36s %10s %10s %10s %10s %10s %10s\n", "counter", "std_us", "std_pp", "std_mm", "lep_us", "lep_pp", "lep_mm");
    for (auto& kv : all) {
        fprintf(ft, "%-36s", kv.first.c_str());
        for (double x : kv.second) fprintf(ft, " %10.0f", x);
        fprintf(ft, "\n");
    }
}

void Analyse() {
    TStopwatch sw;
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(0);
    TH1::AddDirectory(false);
    for (int k = 0; k < 6; ++k) {
        gF[k] = TFile::Open(Form(kPassFmt, k), "READ");
        if (!gF[k] || gF[k]->IsZombie()) { printf("missing pass file %d: run stage 1 first\n", k); return; }
    }
    const std::string td = simeng::kTabDir;
    auto open = [&](const char* n) { return fopen((td + "/" + n).c_str(), "w"); };
    FILE* f1 = open("closure_ratios.txt");
    FILE* f1b = open("closure_ratio_bins.txt");
    PartClosure(f1, f1b);
    fclose(f1); fclose(f1b);
    FILE* f2 = open("closure_partition.txt");
    PartPartition(f2);
    fclose(f2);
    FILE* f3 = open("closure_flags.txt");
    PartFlags(f3);
    fclose(f3);
    FILE* f4 = open("closure_ghost.txt");
    PartGhost(f4);
    fclose(f4);
    FILE* f5 = open("closure_codes.txt");
    PartCodes(f5);
    fclose(f5);
    FILE* f6 = open("closure_counters.txt");
    PartCounters(f6);
    fclose(f6);
    for (int k = 0; k < 6; ++k) { gF[k]->Close(); delete gF[k]; gF[k] = nullptr; }
    printf("stage 2 done in %.0f s\n", sw.RealTime());
}

}  // namespace scl

#endif  // SIM_CLOSURE_ANALYSIS_H
