// mass_ee_ratio_const_fit_exp_vs_sim_channel.C — Channel-filtered version of
// mass_ee_ratio_const_fit_exp_vs_sim.C. Same EXP side (CB-subtracted) but
// SIM side is now restricted to a single physics channel by filtering on
// `sim_geninfo1` (read DIRECTLY from the SMASH source files because the
// existing output_epem_sim.root does not carry the geninfo branches).
//
// Two channels are processed in a single run, each producing its own
// 2-panel canvas (correction factor + trigger efficiency):
//   • sim_geninfo1 == 7051   → π⁰ Dalitz   (π⁰ → e⁺e⁻γ)
//   • sim_geninfo1 == 17051  → η  Dalitz   (η  → e⁺e⁻γ)
//
// EXP side is identical to the reference macro:
//   exp = N_PT3 / (63 · N_PT2)  with CB subtraction via 2·√(N_++ · N_--).
// SIM side reads SMASH source files via TChain, applies the same purity +
// same-vertex gate used by main.cc (sim/main.cc:55-63), the same trigger
// bit extraction (pt3 = bit 13, pt2 = bit 12, not mutually exclusive), and
// additionally requires sim_geninfo1 == channel-value.
//
// Arguments:
//   1) oa_cut_deg   — OA threshold; default 0 ⇒ no cut (applied to BOTH exp & sim)
//   2) use_weights  — true → weighted χ² fit; false → unweighted; default true
//   3) fit_xmin     — fit range lower bound [GeV/c²]; -1 ⇒ full
//   4) fit_xmax     — fit range upper bound [GeV/c²]; -1 ⇒ full
//   5) n_files      — number of SMASH source files to chain; default 3 (max 10)
//
// Usage:
//   root -l -b -q plots/mass_ee_ratio_const_fit_exp_vs_sim_channel.C
//   root -l -b -q 'plots/mass_ee_ratio_const_fit_exp_vs_sim_channel.C(0,true,0.1,0.8,10)'
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TF1.h>
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
#include <limits>
#include <string>
#include <vector>

namespace {

// -------- Binning (same as reference exp-vs-sim macro) --------
constexpr int    kNb   = 20;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.4;

constexpr double kMeMeV = 0.51099895;
constexpr double kD2R   = 1.74532925199432955e-02;

constexpr Color_t kColExp    = kBlue + 1;
constexpr Color_t kColExpFit = kBlue + 2;
constexpr Color_t kColSim    = kRed  + 1;
constexpr Color_t kColSimFit = kRed  + 2;

// -------- Reconstructed m_ee from spherical (p, θ_deg, φ_deg) per leg ----
double m_ee_rec_GeV(double p1, double th1_deg, double ph1_deg,
                    double p2, double th2_deg, double ph2_deg) {
    const double th1 = th1_deg * kD2R, ph1 = ph1_deg * kD2R;
    const double th2 = th2_deg * kD2R, ph2 = ph2_deg * kD2R;
    const double px1 = p1 * std::sin(th1) * std::cos(ph1);
    const double py1 = p1 * std::sin(th1) * std::sin(ph1);
    const double pz1 = p1 * std::cos(th1);
    const double px2 = p2 * std::sin(th2) * std::cos(ph2);
    const double py2 = p2 * std::sin(th2) * std::sin(ph2);
    const double pz2 = p2 * std::cos(th2);
    const double E1 = std::sqrt(p1 * p1 + kMeMeV * kMeMeV);
    const double E2 = std::sqrt(p2 * p2 + kMeMeV * kMeMeV);
    const double dot = px1 * px2 + py1 * py2 + pz1 * pz2;
    const double m2 = 2.0 * kMeMeV * kMeMeV + 2.0 * (E1 * E2 - dot);
    if (m2 < 0.0) return 0.0;
    return std::sqrt(m2) * 1e-3;
}

double openingAngleDeg(double th1_deg, double ph1_deg,
                       double th2_deg, double ph2_deg) {
    const double th1 = th1_deg * kD2R, ph1 = ph1_deg * kD2R;
    const double th2 = th2_deg * kD2R, ph2 = ph2_deg * kD2R;
    const double cosA = std::sin(th1) * std::sin(th2) *
                            std::cos(ph1 - ph2) +
                        std::cos(th1) * std::cos(th2);
    return std::acos(std::max(-1.0, std::min(1.0, cosA))) / kD2R;
}

// -------- Read SMASH list (strip quotes / commas / whitespace) -------
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

// -------- Common helpers (mirrors exp-vs-sim macro) ------------------
TH1D* drawFromNt(TTree* t, const std::string& cut_or_weight,
                 const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kNb, kXmin, kXmax);
    h->Sumw2();
    t->Draw(("m_ee>>" + name).c_str(), cut_or_weight.c_str(), "goff");
    h->SetDirectory(nullptr);
    return h;
}

TH1D* makeCB(TH1D* pp, TH1D* mm, const std::string& name) {
    TH1D* cb = static_cast<TH1D*>(pp->Clone(name.c_str()));
    cb->SetDirectory(nullptr);
    cb->Reset();
    const int nb = pp->GetNbinsX();
    for (int b = 1; b <= nb; ++b) {
        const double np = pp->GetBinContent(b);
        const double ep = pp->GetBinError(b);
        const double nm = mm->GetBinContent(b);
        const double em = mm->GetBinError(b);
        if (np > 0 && nm > 0) {
            const double cval = 2.0 * std::sqrt(np * nm);
            const double rel  = std::sqrt(std::pow(ep / np, 2) +
                                           std::pow(em / nm, 2));
            cb->SetBinContent(b, cval);
            cb->SetBinError  (b, cval * 0.5 * rel);
        }
    }
    return cb;
}

TH1D* makeSig(TH1D* all, TH1D* cb, const std::string& name) {
    TH1D* sig = static_cast<TH1D*>(all->Clone(name.c_str()));
    sig->SetDirectory(nullptr);
    sig->Add(cb, -1.0);
    return sig;
}

TH1D* makeRatio(TH1D* num, TH1D* den, const std::string& name, double scale) {
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

void styleDot(TH1D* h, Color_t color, Style_t marker) {
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.9);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.25);
}

struct FitResult { double a; double e; double chi2_ndf; int ndf; };

FitResult fitConst(TH1D* h, double fit_xmin, double fit_xmax,
                   bool use_weights, Color_t color) {
    const std::string fname = std::string(h->GetName()) + "_fconst";
    TF1* f = new TF1(fname.c_str(), "[0]", fit_xmin, fit_xmax);
    f->SetParName(0, "a");
    f->SetLineColor(color);
    f->SetLineWidth(2);
    f->SetLineStyle(1);
    const std::string opt = std::string("RQ") + (use_weights ? "" : "W");
    h->Fit(f, opt.c_str());
    FitResult r;
    r.a        = f->GetParameter(0);
    r.e        = f->GetParError(0);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
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

// Channel definition for the run.
struct Channel {
    int         geninfo1;
    std::string short_id;          // for file suffix
    std::string title;             // for canvas title and legend
};

}  // anonymous namespace

void mass_ee_ratio_const_fit_exp_vs_sim_channel(double oa_cut_deg  = 0.0,
                                                bool   use_weights = true,
                                                double fit_xmin    = -1.0,
                                                double fit_xmax    = -1.0,
                                                int    n_files     = 3) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);

    const std::vector<Channel> channels = {
        {7051,  "pi0_dalitz", "#pi^{0} Dalitz (geninfo1 = 7051)"},
        {17051, "eta_dalitz", "#eta Dalitz (geninfo1 = 17051)"},
    };

    const bool   custom_range =
        (fit_xmin >= 0.0 && fit_xmin < kXmax) ||
        (fit_xmax >  0.0 && fit_xmax <= kXmax);
    const double fxlo = (fit_xmin >= 0.0 && fit_xmin < kXmax) ? fit_xmin : kXmin;
    const double fxhi = (fit_xmax >  0.0 && fit_xmax <= kXmax) ? fit_xmax : kXmax;

    // ============================================================================
    // EXP — same as reference macro (CB-subtracted PT3/PT2, factor 63)
    // ============================================================================
    TFile* fe_em = TFile::Open("output_epem_exp.root", "READ");
    TFile* fe_pp = TFile::Open("output_epep_exp.root", "READ");
    TFile* fe_mm = TFile::Open("output_emem_exp.root", "READ");
    if (!fe_em || !fe_pp || !fe_mm ||
        fe_em->IsZombie() || fe_pp->IsZombie() || fe_mm->IsZombie()) {
        std::cerr << "Cannot open one of the EXP files\n";
        return;
    }
    TTree* nt_em = dynamic_cast<TTree*>(fe_em->Get("dilepton_nt"));
    TTree* nt_pp = dynamic_cast<TTree*>(fe_pp->Get("dilepton_nt"));
    TTree* nt_mm = dynamic_cast<TTree*>(fe_mm->Get("dilepton_nt"));
    if (!nt_em || !nt_pp || !nt_mm) {
        std::cerr << "dilepton_nt missing in one of the EXP files\n";
        return;
    }

    const std::string oa_exp =
        (oa_cut_deg > 0.0)
            ? std::string(" && oa>") + std::to_string(oa_cut_deg)
            : std::string{};
    const std::string cut_pt3_exp = "trigbit==8192" + oa_exp;
    const std::string cut_pt2_exp = "trigbit==4096" + oa_exp;

    TH1D* em_p3 = drawFromNt(nt_em, cut_pt3_exp, "em_p3_e_ch");
    TH1D* pp_p3 = drawFromNt(nt_pp, cut_pt3_exp, "pp_p3_e_ch");
    TH1D* mm_p3 = drawFromNt(nt_mm, cut_pt3_exp, "mm_p3_e_ch");
    TH1D* em_p2 = drawFromNt(nt_em, cut_pt2_exp, "em_p2_e_ch");
    TH1D* pp_p2 = drawFromNt(nt_pp, cut_pt2_exp, "pp_p2_e_ch");
    TH1D* mm_p2 = drawFromNt(nt_mm, cut_pt2_exp, "mm_p2_e_ch");
    TH1D* cb_p3  = makeCB (pp_p3, mm_p3, "cb_p3_e_ch");
    TH1D* cb_p2  = makeCB (pp_p2, mm_p2, "cb_p2_e_ch");
    TH1D* sig_p3 = makeSig(em_p3, cb_p3, "sig_p3_e_ch");
    TH1D* sig_p2 = makeSig(em_p2, cb_p2, "sig_p2_e_ch");

    TH1D* r_corr_exp = makeRatio(sig_p2, sig_p3, "r_corr_exp_ch", 63.0);
    TH1D* r_eff_exp  = makeRatio(sig_p3, sig_p2, "r_eff_exp_ch",  1.0 / 63.0);
    styleDot(r_corr_exp, kColExp, 20);
    styleDot(r_eff_exp,  kColExp, 20);

    // ============================================================================
    // SIM — chain SMASH source files and fill per-channel PT3/PT2 spectra
    // ============================================================================
    const std::string list_path =
        "/home/damian/hdd1/HADES/pp45/sim/smash/smash_lepton.list";
    auto files = readList(list_path, n_files);
    if (files.empty()) {
        std::cerr << "No SMASH source files found\n";
        return;
    }
    std::cout << "Chaining " << files.size() << " SMASH lepton source file(s):\n";
    TChain ch("EpEm_ID");
    for (const auto& p : files) {
        std::cout << "  " << p << "\n";
        ch.Add(p.c_str());
    }
    const Long64_t n_entries = ch.GetEntries();
    std::cout << "Total entries: " << n_entries << "\n";

    TTreeReader r(&ch);
    TTreeReaderValue<float> v_ep_id  (r, "ep_sim_id");
    TTreeReaderValue<float> v_em_id  (r, "em_sim_id");
    TTreeReaderValue<float> v_ep_g1  (r, "ep_sim_geninfo1");
    TTreeReaderValue<float> v_ep_g2  (r, "ep_sim_geninfo2");
    TTreeReaderValue<float> v_em_g2  (r, "em_sim_geninfo2");
    TTreeReaderValue<float> v_ep_w   (r, "ep_sim_genweight");
    TTreeReaderValue<float> v_trigb  (r, "trigbit");
    TTreeReaderValue<float> v_isBest (r, "isBest");
    TTreeReaderValue<float> v_vz     (r, "eVertReco_z");
    TTreeReaderValue<float> v_ep_p   (r, "ep_p");
    TTreeReaderValue<float> v_ep_th  (r, "ep_theta");
    TTreeReaderValue<float> v_ep_ph  (r, "ep_phi");
    TTreeReaderValue<float> v_em_p   (r, "em_p");
    TTreeReaderValue<float> v_em_th  (r, "em_theta");
    TTreeReaderValue<float> v_em_ph  (r, "em_phi");

    const int n_ch = static_cast<int>(channels.size());
    std::vector<TH1D*> h_p3_sim(n_ch, nullptr);
    std::vector<TH1D*> h_p2_sim(n_ch, nullptr);
    for (int i = 0; i < n_ch; ++i) {
        h_p3_sim[i] = new TH1D(Form("h_p3_sim_g%d", channels[i].geninfo1),
                               "", kNb, kXmin, kXmax);
        h_p3_sim[i]->Sumw2(); h_p3_sim[i]->SetDirectory(nullptr);
        h_p2_sim[i] = new TH1D(Form("h_p2_sim_g%d", channels[i].geninfo1),
                               "", kNb, kXmin, kXmax);
        h_p2_sim[i]->Sumw2(); h_p2_sim[i]->SetDirectory(nullptr);
    }

    std::vector<Long64_t> n_pass(n_ch, 0), n_pt3(n_ch, 0), n_pt2(n_ch, 0);
    Long64_t n_total = 0, n_purity = 0;

    while (r.Next()) {
        ++n_total;
        // Quality + purity gates (mirror main.cc:47-63)
        if (static_cast<int>(*v_isBest) != 1)              continue;
        if (*v_vz < -500.0)                                continue;
        if (static_cast<int>(*v_ep_id) != 2)               continue;
        if (static_cast<int>(*v_em_id) != 3)               continue;
        if (*v_ep_g2 != *v_em_g2)                          continue;
        ++n_purity;

        // OA cut (applied to SIM identically to EXP)
        if (oa_cut_deg > 0.0) {
            const double oa = openingAngleDeg(*v_ep_th, *v_ep_ph,
                                              *v_em_th, *v_em_ph);
            if (oa <= oa_cut_deg) continue;
        }

        const int    g1     = static_cast<int>(*v_ep_g1);
        const int    trigb  = static_cast<int>(*v_trigb);
        const bool   is_pt3 = (trigb >> 13) & 1;
        const bool   is_pt2 = (trigb >> 12) & 1;
        if (!is_pt3 && !is_pt2) continue;            // not interesting

        const double w = use_weights ? static_cast<double>(*v_ep_w) : 1.0;
        const double m_ee = m_ee_rec_GeV(*v_ep_p, *v_ep_th, *v_ep_ph,
                                         *v_em_p, *v_em_th, *v_em_ph);

        for (int i = 0; i < n_ch; ++i) {
            if (g1 != channels[i].geninfo1) continue;
            ++n_pass[i];
            if (is_pt3) { h_p3_sim[i]->Fill(m_ee, w); ++n_pt3[i]; }
            if (is_pt2) { h_p2_sim[i]->Fill(m_ee, w); ++n_pt2[i]; }
        }
    }
    std::cout << "Scanned " << n_total << "; passed purity+vertex "
              << n_purity << " ("
              << (100.0 * n_purity / std::max<Long64_t>(1, n_total))
              << " %)\n";
    for (int i = 0; i < n_ch; ++i) {
        std::cout << "  " << channels[i].title << ":  matched "
                  << n_pass[i] << "   pt3 " << n_pt3[i] << "   pt2 "
                  << n_pt2[i] << "\n";
    }
    std::cout << "Fit mode: " << (use_weights ? "WEIGHTED" : "UNWEIGHTED")
              << "   Fit range: [" << fxlo << ", " << fxhi << "] GeV/c^2"
              << (custom_range ? "  (custom)" : "  (full)") << "\n\n";

    // EXP fits (shared across channels — only sim differs)
    std::cout << "EXP constant fits (all channels):\n";
    const FitResult f_corr_exp =
        fitConst(r_corr_exp, fxlo, fxhi, use_weights, kColExpFit);
    const FitResult f_eff_exp =
        fitConst(r_eff_exp,  fxlo, fxhi, use_weights, kColExpFit);
    std::cout << "  corr  a = " << f_corr_exp.a << " ± " << f_corr_exp.e
              << "  chi2/ndf = " << (f_corr_exp.chi2_ndf * f_corr_exp.ndf)
              << "/" << f_corr_exp.ndf << "\n"
              << "  eff   a = " << f_eff_exp.a << " ± " << f_eff_exp.e
              << "  chi2/ndf = " << (f_eff_exp.chi2_ndf * f_eff_exp.ndf)
              << "/" << f_eff_exp.ndf << "\n";

    // ============================================================================
    // Per-channel canvases (correction factor + trigger efficiency)
    // ============================================================================
    gSystem->mkdir("plots/output", true);
    const std::string oa_suffix =
        (oa_cut_deg > 0.0) ? Form("_oa%g", oa_cut_deg) : std::string{};
    const std::string wt_suffix = use_weights ? std::string{} : "_unw";
    const std::string fr_suffix =
        custom_range ? Form("_fit%g-%g", fxlo, fxhi) : std::string{};
    const std::string nf_suffix = Form("_n%d", (int)files.size());

    for (int i = 0; i < n_ch; ++i) {
        TH1D* r_corr_sim = makeRatio(h_p2_sim[i], h_p3_sim[i],
                                     Form("r_corr_sim_g%d", channels[i].geninfo1),
                                     1.0);
        TH1D* r_eff_sim  = makeRatio(h_p3_sim[i], h_p2_sim[i],
                                     Form("r_eff_sim_g%d", channels[i].geninfo1),
                                     1.0);
        styleDot(r_corr_sim, kColSim, 21);
        styleDot(r_eff_sim,  kColSim, 21);

        const std::string oa_tag =
            (oa_cut_deg > 0.0)
                ? Form(", OA > %g#circ", oa_cut_deg)
                : std::string{};

        // Reset titles per channel — preserves the same axis labels as ref macro.
        r_corr_exp->SetTitle((std::string(
            "Trigger correction factor — ") + channels[i].title +
            "  (exp vs sim, 20 fixed bins" + oa_tag + ");"
            "M_{e^{+}e^{-}} [GeV/c^{2}];"
            "63 #upoint N_{PT2}/N_{PT3} (exp)   N_{PT2}/N_{PT3} (sim)").c_str());
        r_eff_exp->SetTitle((std::string(
            "Trigger efficiency — ") + channels[i].title +
            "  (exp vs sim, 20 fixed bins" + oa_tag + ");"
            "M_{e^{+}e^{-}} [GeV/c^{2}];"
            "N_{PT3}/(63#upointN_{PT2}) (exp)   N_{PT3}/N_{PT2} (sim)").c_str());

        const FitResult f_corr_sim =
            fitConst(r_corr_sim, fxlo, fxhi, use_weights, kColSimFit);
        const FitResult f_eff_sim =
            fitConst(r_eff_sim,  fxlo, fxhi, use_weights, kColSimFit);
        std::cout << "[" << channels[i].title << "] SIM fits:\n"
                  << "  corr  a = " << f_corr_sim.a << " ± " << f_corr_sim.e
                  << "  chi2/ndf = " << (f_corr_sim.chi2_ndf * f_corr_sim.ndf)
                  << "/" << f_corr_sim.ndf << "\n"
                  << "  eff   a = " << f_eff_sim.a << " ± " << f_eff_sim.e
                  << "  chi2/ndf = " << (f_eff_sim.chi2_ndf * f_eff_sim.ndf)
                  << "/" << f_eff_sim.ndf << "\n";

        TCanvas* c = new TCanvas(
            Form("c_chan_%s", channels[i].short_id.c_str()),
            Form("Exp vs Sim trigger ratios — %s",
                 channels[i].title.c_str()),
            1500, 600);
        c->Divide(2, 1, 0.001, 0.001);

        auto drawComparePad = [&](int idx, TH1D* h_exp, TH1D* h_sim,
                                  const FitResult& f_exp,
                                  const FitResult& f_sim,
                                  double ref_y, double y_max_floor) {
            c->cd(idx);
            const double pad_lm = 0.13, pad_rm = 0.04, pad_bm = 0.13, pad_tm = 0.10;
            gPad->SetMargin(pad_lm, pad_rm, pad_bm, pad_tm);
            const double ymax = std::max(safeMaxAbs(h_exp), safeMaxAbs(h_sim));
            const double y_hi = std::max(ymax * 1.30, y_max_floor);
            h_exp->GetYaxis()->SetRangeUser(0.0, y_hi);
            h_exp->Draw("E1");
            h_sim->Draw("E1 SAME");

            TLine* lref = new TLine(kXmin, ref_y, kXmax, ref_y);
            lref->SetLineStyle(3); lref->SetLineColor(kGray + 2); lref->Draw();

            TLatex tex;
            tex.SetNDC();
            tex.SetTextAlign(11);
            tex.SetTextSize(0.040);
            tex.SetTextColor(kColExpFit);
            tex.DrawLatex(0.16, 0.85,
                Form("exp:  a = %.4f #pm %.4f", f_exp.a, f_exp.e));
            tex.SetTextSize(0.030);
            tex.DrawLatex(0.16, 0.81,
                Form("       #chi^{2}/ndf = %.1f / %d",
                     f_exp.chi2_ndf * f_exp.ndf, f_exp.ndf));

            tex.SetTextAlign(31);
            tex.SetTextSize(0.040);
            tex.SetTextColor(kColSimFit);
            tex.DrawLatex(0.95, 0.85,
                Form("sim:  a = %.4f #pm %.4f", f_sim.a, f_sim.e));
            tex.SetTextSize(0.030);
            tex.DrawLatex(0.95, 0.81,
                Form("#chi^{2}/ndf = %.1f / %d       ",
                     f_sim.chi2_ndf * f_sim.ndf, f_sim.ndf));

            TLegend* leg = new TLegend(0.40, 0.74, 0.65, 0.83);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
            leg->AddEntry(h_exp, "exp",  "lpe");
            leg->AddEntry(h_sim, "sim",  "lpe");
            leg->Draw();
        };

        drawComparePad(1, r_corr_exp, r_corr_sim, f_corr_exp, f_corr_sim,
                       1.0, 3.5);
        drawComparePad(2, r_eff_exp,  r_eff_sim,  f_eff_exp,  f_eff_sim,
                       1.0, 1.2);

        const std::string base =
            "plots/output/mass_ee_ratio_const_fit_exp_vs_sim_" +
            channels[i].short_id +
            oa_suffix + wt_suffix + fr_suffix + nf_suffix;
        c->SaveAs((base + ".pdf").c_str());
        c->SaveAs((base + ".png").c_str());
        std::cout << "Saved: " << base << ".{pdf,png}\n";
    }

    fe_em->Close();
    fe_pp->Close();
    fe_mm->Close();
}
