// mass_ee_pt_ratio_final.C
// =========================================================================
// EXP analysis: 63 * N_PT2 / N_PT3 CB-subtracted signal ratio with a
// four-region variable binning layout (empirical for smoothness):
//
//   [0.00, 0.05] : 20 bins x 2.5  MeV   (DENSE LOW, pi0 region)
//   [0.05, 0.14] : 3  bins x 30   MeV   (LOW-FADE, Config U3)
//   [0.14, 0.80] : 33 bins x 20   MeV   (MID, uniform 20 MeV)
//   [0.80, 1.40] : 2  bins variable     (HIGH, Config F: 0.80,0.95,1.40)
//   TOTAL: 58 bins
//
// HIGH-region empirical choice:
//   Edges {0.80, 0.95, 1.40} (Config F) is the ONLY tested layout in
//   that region that produces no NEGATIVE/zero ratio bins — the
//   high-mass tail above ~1.10 has CB > epem so PT2 signal goes
//   negative there; the wide 0.95-1.40 bin absorbs it cleanly.
//
// LOW-FADE empirical choice:
//   Three 30 MeV bins on [0.05, 0.14] (Config U3) is the chosen
//   layout from the _explore_low_fade_005 scan. It produces zero BAD,
//   monotonic descent (1.929 -> 1.870 -> 1.634), and the lowest mean
//   per-bin sigma (0.043) among configs with >=3 bins. The single
//   "JUMP" at the last bin reflects the genuine physical transition
//   into the MID pol0 plateau (~1.82) and is expected, not noise.
//
// pol0 fits:
//   C_MID = pol0 over MID-only range [0.14, 0.80]
//   C_EXT = pol0 over EXTENDED range [0.05, 0.80] (includes LOW-FADE)
//   Reports both, plus DeltaC and significance.
//
// Trigger cuts ONLY (everything else baked at main.cc EXP level):
//   PT3 = trigbit==8192
//   PT2 = trigbit==4096
//
// Output:
//   plots/output/mass_ee_pt_ratio_final.{pdf,png}
//
// Usage: root -l -b -q plots/mass_ee_pt_ratio_final.C
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

namespace {

constexpr const char* kFileEm = "output_epem_exp.root";
constexpr const char* kFilePp = "output_epep_exp.root";
constexpr const char* kFileMm = "output_emem_exp.root";
constexpr const char* kTree   = "dilepton_nt";

constexpr int    kNb   = 58;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.4;

// Trigger downscale: PT2 is downscaled by 64, so we multiply by 63.
constexpr double kTrigCorr = 63.0;

// Trigger bit cuts (EXP) — NO additional selection.
constexpr const char* kCutPT3 = "trigbit==8192";
constexpr const char* kCutPT2 = "trigbit==4096";

// MID-only pol0 fit range.
constexpr double kMdLo = 0.14;
constexpr double kMdHi = 0.80;

// EXTENDED pol0 fit range (includes LOW-FADE region).
constexpr double kExtLo = 0.05;
constexpr double kExtHi = 0.80;

// Four-region transition points.
constexpr double kDenseHi = 0.05;  // dense low ends here
constexpr double kLowHi   = 0.14;  // low-fade ends, MID starts
constexpr double kHighLo  = 0.80;  // MID ends, HIGH starts

// Build the variable-width bin-edge array (58 bins total).
//   [0.00, 0.05] : 20 bins x 2.5  MeV   (DENSE LOW)
//   [0.05, 0.14] : 3  bins x 30   MeV   (LOW-FADE, Config U3)
//   [0.14, 0.80] : 33 bins x 20   MeV   (MID)
//   [0.80, 1.40] : 2  bins variable     (HIGH, Config F: 0.95 split)
const double* getEdges() {
    static double edges[kNb + 1];
    static bool   built = false;
    if (!built) {
        int k = 0;
        // DENSE LOW: 0.00, 0.0025, ..., 0.05 (21 edges => 20 bins)
        for (int i = 0; i <= 20; ++i) edges[k++] = i * 0.0025;
        // LOW-FADE (U3): 3 bins x 30 MeV up to 0.14 (3 edges added)
        edges[k++] = 0.0800;
        edges[k++] = 0.1100;
        edges[k++] = 0.1400;
        // MID: 33 bins x 20 MeV up to 0.80 (33 edges added)
        for (int i = 1; i <= 33; ++i) edges[k++] = 0.140 + i * 0.020;
        // HIGH Config F: 0.95, 1.40 (2 edges added => 2 bins)
        edges[k++] = 0.950;
        edges[k++] = 1.400;
        // sanity: edges[k-1] should be 1.400 and k should be kNb+1 = 59
        built = true;
    }
    return edges;
}

// Fill a histogram from a TTree using TTree::Draw with cut.
TH1D* drawWithCut(TTree* t, const char* expr,
                  const char* cut, const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kNb, getEdges());
    h->Sumw2();
    t->Draw((std::string(expr) + ">>" + name).c_str(), cut, "goff");
    h->SetDirectory(nullptr);
    return h;
}

// CB = 2*sqrt(N_++ * N_--) per bin with proper error propagation.
TH1D* makeCB(TH1D* h_pp, TH1D* h_mm, const std::string& name) {
    TH1D* cb = static_cast<TH1D*>(h_pp->Clone(name.c_str()));
    cb->SetDirectory(nullptr);
    cb->Reset();
    for (int b = 1; b <= h_pp->GetNbinsX(); ++b) {
        const double n_pp = h_pp->GetBinContent(b);
        const double e_pp = h_pp->GetBinError(b);
        const double n_mm = h_mm->GetBinContent(b);
        const double e_mm = h_mm->GetBinError(b);
        if (n_pp <= 0.0 || n_mm <= 0.0 ||
            !std::isfinite(n_pp) || !std::isfinite(n_mm)) {
            cb->SetBinContent(b, 0.0);
            cb->SetBinError  (b, 0.0);
            continue;
        }
        const double val = 2.0 * std::sqrt(n_pp * n_mm);
        const double rel_pp = e_pp / n_pp;
        const double rel_mm = e_mm / n_mm;
        const double err    = val * 0.5 *
            std::sqrt(rel_pp * rel_pp + rel_mm * rel_mm);
        cb->SetBinContent(b, val);
        cb->SetBinError  (b, err);
    }
    return cb;
}

// Signal = epem - CB with quadrature errors.
TH1D* makeSignal(TH1D* h_em, TH1D* h_cb, const std::string& name) {
    TH1D* sig = static_cast<TH1D*>(h_em->Clone(name.c_str()));
    sig->SetDirectory(nullptr);
    sig->Reset();
    for (int b = 1; b <= h_em->GetNbinsX(); ++b) {
        const double n_em = h_em->GetBinContent(b);
        const double e_em = h_em->GetBinError(b);
        const double n_cb = h_cb->GetBinContent(b);
        const double e_cb = h_cb->GetBinError(b);
        const double val  = n_em - n_cb;
        const double err  = std::sqrt(e_em * e_em + e_cb * e_cb);
        sig->SetBinContent(b, val);
        sig->SetBinError  (b, err);
    }
    return sig;
}

// Ratio = factor * num / den per bin with full error propagation.
// Skip bins with den <= 0 (set content & error to 0).
TH1D* makeRatio(TH1D* num, TH1D* den, double factor,
                const std::string& name) {
    TH1D* r = static_cast<TH1D*>(num->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    for (int b = 1; b <= num->GetNbinsX(); ++b) {
        const double n  = num->GetBinContent(b);
        const double en = num->GetBinError(b);
        const double d  = den->GetBinContent(b);
        const double ed = den->GetBinError(b);
        if (d <= 0.0 || !std::isfinite(d) || !std::isfinite(n)) {
            r->SetBinContent(b, 0.0);
            r->SetBinError  (b, 0.0);
            continue;
        }
        const double val = factor * n / d;
        const double rel = std::sqrt(
            (n != 0.0 ? std::pow(en / n, 2) : 0.0) +
            std::pow(ed / d, 2));
        r->SetBinContent(b, val);
        r->SetBinError  (b, std::abs(val) * rel);
    }
    return r;
}

// Width-weighted average ratio over a region [xlo, xhi]:
//   mean = sum(content_i * width_i) / sum(width_i)
// Uncertainty propagated as:
//   sigma = sqrt(sum((error_i * width_i)^2)) / sum(width_i)
// Skips bins with zero content AND zero error (no data).
void regionMean(TH1D* h, double xlo, double xhi,
                double& mean, double& sigma, int& nbins_used) {
    double sum_vw   = 0.0;
    double sum_w    = 0.0;
    double sum_e2w2 = 0.0;
    nbins_used = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double lo = h->GetXaxis()->GetBinLowEdge(b);
        const double hi = h->GetXaxis()->GetBinUpEdge (b);
        if (lo < xlo - 1e-9) continue;
        if (hi > xhi + 1e-9) continue;
        const double w = h->GetBinWidth(b);
        const double v = h->GetBinContent(b);
        const double e = h->GetBinError  (b);
        // Skip bins that carry NO information (content and error both zero).
        if (v == 0.0 && e == 0.0) continue;
        sum_vw   += v * w;
        sum_w    += w;
        sum_e2w2 += (e * w) * (e * w);
        ++nbins_used;
    }
    if (sum_w > 0.0) {
        mean  = sum_vw / sum_w;
        sigma = std::sqrt(sum_e2w2) / sum_w;
    } else {
        mean = 0.0;
        sigma = 0.0;
    }
}

void styleHist(TH1D* h, Color_t c, Style_t m) {
    h->SetLineColor(c);
    h->SetMarkerColor(c);
    h->SetMarkerStyle(m);
    h->SetMarkerSize(0.8);
    h->SetLineWidth(1);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.040);
    h->GetYaxis()->SetLabelSize(0.040);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.15);
}

}  // anonymous namespace

void mass_ee_pt_ratio_innermdc() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0.0);

    // ---------- Open input files ----------
    TFile* f_em = TFile::Open(kFileEm, "READ");
    TFile* f_pp = TFile::Open(kFilePp, "READ");
    TFile* f_mm = TFile::Open(kFileMm, "READ");
    if (!f_em || f_em->IsZombie() ||
        !f_pp || f_pp->IsZombie() ||
        !f_mm || f_mm->IsZombie()) {
        std::cerr << "Cannot open one of the EXP input files.\n";
        return;
    }
    TTree* t_em = dynamic_cast<TTree*>(f_em->Get(kTree));
    TTree* t_pp = dynamic_cast<TTree*>(f_pp->Get(kTree));
    TTree* t_mm = dynamic_cast<TTree*>(f_mm->Get(kTree));
    if (!t_em || !t_pp || !t_mm) {
        std::cerr << "Tree " << kTree << " missing in one of the EXP files.\n";
        return;
    }

    std::cout << "Input:\n";
    std::cout << "  epem : " << kFileEm << " ("
              << t_em->GetEntries() << " entries)\n";
    std::cout << "  epep : " << kFilePp << " ("
              << t_pp->GetEntries() << " entries)\n";
    std::cout << "  emem : " << kFileMm << " ("
              << t_mm->GetEntries() << " entries)\n";
    std::cout << "Binning (four-region empirical, " << kNb << " bins total):\n";
    std::cout << "  [0.00, 0.05] : 20 x 2.5  MeV (DENSE LOW)\n";
    std::cout << "  [0.05, 0.14] : 3  x 30   MeV (LOW-FADE, Config U3)\n";
    std::cout << "  [0.14, 0.80] : 33 x 20   MeV (MID)\n";
    std::cout << "  [0.80, 1.40] : 2 variable    (HIGH, Config F: 0.80/0.95/1.40)\n";
    std::cout << "Trigger cuts: PT3='" << kCutPT3 << "'  PT2='"
              << kCutPT2 << "'\n";
    std::cout << "Correction factor (PT2 downscale 64 -> factor): "
              << kTrigCorr << "\n\n";

    // Sanity check edges. Layout: 20 dense + 3 low-fade + 33 mid + 2 high.
    // Indices (edge index): 0=>0.000, 20=>0.050, 23=>0.140, 56=>0.800,
    // 57=>0.950, 58=>1.400.
    const double* edges = getEdges();
    std::printf("Edge sanity: edges[0]=%.4f  edges[20]=%.4f  "
                "edges[23]=%.4f  edges[56]=%.4f  "
                "edges[57]=%.4f  edges[58]=%.4f\n",
                edges[0], edges[20], edges[23], edges[56],
                edges[57], edges[58]);

    gSystem->mkdir("plots/output", true);

    // ---------- Build raw mass histograms per trigger ----------
    TH1D* h_epem_PT3 = drawWithCut(t_em, "m_ee", kCutPT3, "h_epem_PT3");
    TH1D* h_epem_PT2 = drawWithCut(t_em, "m_ee", kCutPT2, "h_epem_PT2");
    TH1D* h_epep_PT3 = drawWithCut(t_pp, "m_ee", kCutPT3, "h_epep_PT3");
    TH1D* h_epep_PT2 = drawWithCut(t_pp, "m_ee", kCutPT2, "h_epep_PT2");
    TH1D* h_emem_PT3 = drawWithCut(t_mm, "m_ee", kCutPT3, "h_emem_PT3");
    TH1D* h_emem_PT2 = drawWithCut(t_mm, "m_ee", kCutPT2, "h_emem_PT2");

    std::cout << "Raw yields:\n";
    std::cout << "  PT3: epem=" << h_epem_PT3->Integral()
              << "  ++="       << h_epep_PT3->Integral()
              << "  --="       << h_emem_PT3->Integral() << "\n";
    std::cout << "  PT2: epem=" << h_epem_PT2->Integral()
              << "  ++="       << h_epep_PT2->Integral()
              << "  --="       << h_emem_PT2->Integral() << "\n";

    // ---------- CB per trigger ----------
    TH1D* h_CB_PT3 = makeCB(h_epep_PT3, h_emem_PT3, "h_CB_PT3");
    TH1D* h_CB_PT2 = makeCB(h_epep_PT2, h_emem_PT2, "h_CB_PT2");

    // ---------- Signal per trigger ----------
    TH1D* h_sig_PT3 = makeSignal(h_epem_PT3, h_CB_PT3, "h_sig_PT3");
    TH1D* h_sig_PT2 = makeSignal(h_epem_PT2, h_CB_PT2, "h_sig_PT2");

    std::cout << "Signal yields (epem - CB):\n";
    std::cout << "  PT3 signal = " << h_sig_PT3->Integral() << "\n";
    std::cout << "  PT2 signal = " << h_sig_PT2->Integral() << "\n";

    // ---------- Ratio = 63 * sig_PT2 / sig_PT3 ----------
    TH1D* h_ratio = makeRatio(h_sig_PT2, h_sig_PT3, kTrigCorr,
                              "h_ratio_final");
    styleHist(h_ratio, kBlack, 20);

    // ---------- MID-only pol0 fit over [0.14, 0.80] ----------
    TF1* f_mid = new TF1("pol0_mid", "[0]", kMdLo, kMdHi);
    h_ratio->Fit(f_mid, "RQ0");
    double C       = f_mid->GetParameter(0);
    double C_err   = f_mid->GetParError(0);
    double chi2    = f_mid->GetChisquare();
    int    ndf     = f_mid->GetNDF();
    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;

    // Clean any auto-attached fit so we draw it ourselves with control.
    if (h_ratio->GetListOfFunctions())
        h_ratio->GetListOfFunctions()->Clear();

    // ---------- EXTENDED pol0 fit over [0.05, 0.80] (includes LOW-FADE) ----
    TF1* f_ext = new TF1("pol0_ext", "[0]", kExtLo, kExtHi);
    h_ratio->Fit(f_ext, "RQ0");
    double C_ext     = f_ext->GetParameter(0);
    double C_ext_err = f_ext->GetParError(0);
    double chi2_ext  = f_ext->GetChisquare();
    int    ndf_ext   = f_ext->GetNDF();
    double chi2ndf_ext = (ndf_ext > 0) ? chi2_ext / ndf_ext : 0.0;

    // Clean again after second fit.
    if (h_ratio->GetListOfFunctions())
        h_ratio->GetListOfFunctions()->Clear();

    // ---------- UNWEIGHTED pol0 fit over [0.05, 0.80] (every bin w=1) ------
    // "W" option: ignore per-bin errors, treat every bin as equal weight.
    TF1* f_uw = new TF1("pol0_uw", "[0]", kExtLo, kExtHi);
    h_ratio->Fit(f_uw, "RQ0W");
    double C_uw     = f_uw->GetParameter(0);
    double C_uw_err = f_uw->GetParError(0);
    double chi2_uw  = f_uw->GetChisquare();
    int    ndf_uw   = f_uw->GetNDF();
    double chi2ndf_uw = (ndf_uw > 0) ? chi2_uw / ndf_uw : 0.0;

    // Clean again after third fit.
    if (h_ratio->GetListOfFunctions())
        h_ratio->GetListOfFunctions()->Clear();

    // Manual unweighted mean cross-check over bin CENTERS in [kExtLo, kExtHi].
    // Skip bins that carry NO information (content and error both zero).
    double uw_sum = 0.0;
    int    uw_n   = 0;
    for (int b = 1; b <= h_ratio->GetNbinsX(); ++b) {
        const double xc = h_ratio->GetXaxis()->GetBinCenter(b);
        if (xc < kExtLo - 1e-9 || xc > kExtHi + 1e-9) continue;
        const double v = h_ratio->GetBinContent(b);
        const double e = h_ratio->GetBinError  (b);
        if (v == 0.0 && e == 0.0) continue;
        uw_sum += v;
        ++uw_n;
    }
    double uw_mean_manual = (uw_n > 0) ? uw_sum / uw_n : 0.0;
    double uw_var = 0.0;
    for (int b = 1; b <= h_ratio->GetNbinsX(); ++b) {
        const double xc = h_ratio->GetXaxis()->GetBinCenter(b);
        if (xc < kExtLo - 1e-9 || xc > kExtHi + 1e-9) continue;
        const double v = h_ratio->GetBinContent(b);
        const double e = h_ratio->GetBinError  (b);
        if (v == 0.0 && e == 0.0) continue;
        const double d = v - uw_mean_manual;
        uw_var += d * d;
    }
    uw_var = (uw_n > 1) ? uw_var / (uw_n - 1) : 0.0;
    double uw_sigma_manual = (uw_n > 0) ? std::sqrt(uw_var / uw_n) : 0.0;

    // Difference and significance.
    const double dC        = C_ext - C;
    const double dC_sigma  = std::sqrt(C_err * C_err +
                                       C_ext_err * C_ext_err);
    const double dC_signif = (dC_sigma > 0.0) ?
                             std::abs(dC) / dC_sigma : 0.0;

    // ---------- Per-region simple width-weighted means ----------
    double mean_low,  sig_low;   int nb_low  = 0;
    double mean_mid,  sig_mid;   int nb_mid  = 0;
    double mean_high, sig_high;  int nb_high = 0;
    regionMean(h_ratio, 0.00,    kLowHi,  mean_low,  sig_low,  nb_low);
    regionMean(h_ratio, kLowHi,  kHighLo, mean_mid,  sig_mid,  nb_mid);
    regionMean(h_ratio, kHighLo, kXmax,   mean_high, sig_high, nb_high);

    // ---------- Per-bin HIGH-region report ----------
    std::printf("\n========== HIGH REGION PER-BIN RATIOS ==========\n");
    for (int b = 1; b <= h_ratio->GetNbinsX(); ++b) {
        const double lo = h_ratio->GetXaxis()->GetBinLowEdge(b);
        const double hi = h_ratio->GetXaxis()->GetBinUpEdge (b);
        if (lo < kHighLo - 1e-9) continue;
        const double v = h_ratio->GetBinContent(b);
        const double e = h_ratio->GetBinError  (b);
        std::printf("  HIGH bin %d [%.3f, %.3f] (width=%.3f): "
                    "r = %.4f +- %.4f\n",
                    b, lo, hi, hi - lo, v, e);
    }
    std::printf("\n========== LOW-FADE BIN [%.2f, %.2f] ==========\n",
                kDenseHi, kLowHi);
    for (int b = 1; b <= h_ratio->GetNbinsX(); ++b) {
        const double lo = h_ratio->GetXaxis()->GetBinLowEdge(b);
        const double hi = h_ratio->GetXaxis()->GetBinUpEdge (b);
        if (lo < kDenseHi - 1e-9) continue;
        if (hi > kLowHi   + 1e-9) continue;
        const double v = h_ratio->GetBinContent(b);
        const double e = h_ratio->GetBinError  (b);
        std::printf("  LOW-FADE bin %d [%.3f, %.3f] (width=%.3f): "
                    "r = %.4f +- %.4f\n",
                    b, lo, hi, hi - lo, v, e);
    }

    std::printf("\n========== pol0 MID FIT [%.2f, %.2f] ==========\n",
                kMdLo, kMdHi);
    std::printf("  C_MID    = %.4f +- %.4f\n", C, C_err);
    std::printf("  chi2/ndf = %.2f / %d = %.3f\n", chi2, ndf, chi2ndf);

    std::printf("\n========== pol0 EXTENDED FIT [%.2f, %.2f] ==========\n",
                kExtLo, kExtHi);
    std::printf("  C_EXT    = %.4f +- %.4f\n", C_ext, C_ext_err);
    std::printf("  chi2/ndf = %.2f / %d = %.3f\n",
                chi2_ext, ndf_ext, chi2ndf_ext);

    std::printf("\n========== DELTA C ==========\n");
    std::printf("  DeltaC = C_EXT - C_MID = %.4f - %.4f = %.4f\n",
                C_ext, C, dC);
    std::printf("  sigma_DeltaC = sqrt(s_MID^2 + s_EXT^2) = %.4f\n",
                dC_sigma);
    std::printf("  |DeltaC| / sigma = %.3f sigma\n", dC_signif);

    // ---------- UNWEIGHTED report ----------
    std::printf("\n========== pol0 UNWEIGHTED FIT [%.2f, %.2f] ==========\n",
                kExtLo, kExtHi);
    std::printf("  C_UW     = %.4f +- %.4f  (TH1::Fit with 'W')\n",
                C_uw, C_uw_err);
    std::printf("  chi2/ndf = %.2f / %d = %.3f\n",
                chi2_uw, ndf_uw, chi2ndf_uw);
    std::printf("  Manual cross-check: <r>_unw = %.4f +- %.4f over %d bins\n",
                uw_mean_manual, uw_sigma_manual, uw_n);

    // Deltas vs the two weighted fits.
    const double dC_uw_mid       = C_uw - C;
    const double dC_uw_mid_sigma = std::sqrt(C_uw_err * C_uw_err +
                                             C_err    * C_err);
    const double dC_uw_mid_nsig  = (dC_uw_mid_sigma > 0.0) ?
                                   std::abs(dC_uw_mid) / dC_uw_mid_sigma : 0.0;

    const double dC_uw_ext       = C_uw - C_ext;
    const double dC_uw_ext_sigma = std::sqrt(C_uw_err  * C_uw_err +
                                             C_ext_err * C_ext_err);
    const double dC_uw_ext_nsig  = (dC_uw_ext_sigma > 0.0) ?
                                   std::abs(dC_uw_ext) / dC_uw_ext_sigma : 0.0;

    std::printf("\n========== DELTA C_UW vs others ==========\n");
    std::printf("  C_UW - C_MID = %+.4f  sigma_comb = %.4f  -> %.2f sigma\n",
                dC_uw_mid, dC_uw_mid_sigma, dC_uw_mid_nsig);
    std::printf("  C_UW - C_EXT = %+.4f  sigma_comb = %.4f  -> %.2f sigma\n",
                dC_uw_ext, dC_uw_ext_sigma, dC_uw_ext_nsig);

    std::printf("\n========== PER-REGION WIDTH-WEIGHTED MEANS ==========\n");
    std::printf("  LOW  [0.00, %.2f]: <ratio> = %.4f +- %.4f  (%d bins used)\n",
                kLowHi, mean_low, sig_low, nb_low);
    std::printf("  MID  [%.2f, %.2f]: <ratio> = %.4f +- %.4f  (%d bins used)\n",
                kLowHi, kHighLo, mean_mid, sig_mid, nb_mid);
    std::printf("  HIGH [%.2f, %.2f]: <ratio> = %.4f +- %.4f  (%d bins used)\n",
                kHighLo, kXmax, mean_high, sig_high, nb_high);

    std::printf("\nComparison vs v2_flatness_subwindow_test:\n");
    std::printf("  previous C_global = 1.835 +- 0.038 on [0.15, 0.80]\n");
    std::printf("  this     C_MID    = %.4f +- %.4f on [%.2f, %.2f]\n",
                C, C_err, kMdLo, kMdHi);
    std::printf("  this     C_EXT    = %.4f +- %.4f on [%.2f, %.2f]\n",
                C_ext, C_ext_err, kExtLo, kExtHi);

    // ---------- Canvas ----------
    TCanvas* c = new TCanvas("c_mass_ee_pt_ratio_innermdc",
                             "63 * N_PT2 / N_PT3 final ratio",
                             1100, 700);
    c->SetMargin(0.11, 0.04, 0.13, 0.09);

    h_ratio->SetTitle(
        "63#upointN_{PT2}/N_{PT3} signal ratio (CB-subtracted), "
        "variable-width HIGH (empirical for smoothness);"
        "M_{e^{+}e^{-}} [GeV/c^{2}];"
        "63 #upoint N_{PT2}^{sig} / N_{PT3}^{sig}");
    h_ratio->GetXaxis()->SetRangeUser(kXmin, kXmax);
    h_ratio->GetYaxis()->SetRangeUser(0.5, 4.0);

    h_ratio->Draw("E1");

    // Horizontal dashed gray reference line at y = 1.
    TLine* l_ref = new TLine(kXmin, 1.0, kXmax, 1.0);
    l_ref->SetLineStyle(2);
    l_ref->SetLineColor(kGray + 2);
    l_ref->SetLineWidth(1);
    l_ref->Draw();

    // Vertical dotted gray lines at binning/fit transitions:
    //   - dense->low-fade at 0.05
    //   - low-fade->mid at 0.14
    //   - mid->high at 0.80
    //   - HIGH internal split at 0.95
    const double trans[] = {kDenseHi, kLowHi, kHighLo, 0.95};
    for (double x : trans) {
        TLine* lt = new TLine(x, 0.5, x, 4.0);
        lt->SetLineStyle(3);
        lt->SetLineColor(kGray + 1);
        lt->SetLineWidth(1);
        lt->Draw();
    }

    // Translucent BLUE error band around C_EXT in [kExtLo, kExtHi].
    // Drawn FIRST so the narrower red MID band sits on top of it.
    TBox* band_ext = new TBox(kExtLo, C_ext - C_ext_err,
                              kExtHi, C_ext + C_ext_err);
    band_ext->SetFillColorAlpha(kBlue, 0.18);
    band_ext->SetLineColor(kBlue);
    band_ext->SetLineStyle(0);
    band_ext->SetFillStyle(1001);
    band_ext->Draw("same");

    // Translucent GREEN error band around C_UW in [kExtLo, kExtHi].
    TBox* band_uw = new TBox(kExtLo, C_uw - C_uw_err,
                             kExtHi, C_uw + C_uw_err);
    band_uw->SetFillColorAlpha(kGreen + 2, 0.22);
    band_uw->SetLineColor(kGreen + 2);
    band_uw->SetLineStyle(0);
    band_uw->SetFillStyle(1001);
    band_uw->Draw("same");

    // Translucent RED error band around C_MID in [kMdLo, kMdHi].
    TBox* band_mid = new TBox(kMdLo, C - C_err, kMdHi, C + C_err);
    band_mid->SetFillColorAlpha(kRed, 0.28);
    band_mid->SetLineColor(kRed);
    band_mid->SetLineStyle(0);
    band_mid->SetFillStyle(1001);
    band_mid->Draw("same");

    // Solid BLUE line for C_EXT over [kExtLo, kExtHi] (dashed style).
    TLine* l_fit_ext = new TLine(kExtLo, C_ext, kExtHi, C_ext);
    l_fit_ext->SetLineColor(kBlue);
    l_fit_ext->SetLineWidth(2);
    l_fit_ext->SetLineStyle(2);
    l_fit_ext->Draw();

    // Solid GREEN line for C_UW over [kExtLo, kExtHi].
    TLine* l_fit_uw = new TLine(kExtLo, C_uw, kExtHi, C_uw);
    l_fit_uw->SetLineColor(kGreen + 2);
    l_fit_uw->SetLineWidth(2);
    l_fit_uw->SetLineStyle(1);
    l_fit_uw->Draw();

    // Solid RED line for C_MID over [kMdLo, kMdHi].
    TLine* l_fit_mid = new TLine(kMdLo, C, kMdHi, C);
    l_fit_mid->SetLineColor(kRed);
    l_fit_mid->SetLineWidth(2);
    l_fit_mid->SetLineStyle(1);
    l_fit_mid->Draw();

    // Re-draw markers on top so points/error bars sit above bands/lines.
    h_ratio->Draw("E1 SAME");

    // Legend.
    TLegend* leg = new TLegend(0.46, 0.66, 0.96, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.026);
    leg->AddEntry(h_ratio, "data (CB-subtracted)", "lpe");
    leg->AddEntry(l_fit_mid,
        Form("pol0 MID [%.2f, %.2f] weighted: C_{MID} = %.3f #pm %.3f, "
             "#chi^{2}/ndf = %.2f",
             kMdLo, kMdHi, C, C_err, chi2ndf), "l");
    leg->AddEntry(l_fit_ext,
        Form("pol0 EXT [%.2f, %.2f] weighted: C_{EXT} = %.3f #pm %.3f, "
             "#chi^{2}/ndf = %.2f",
             kExtLo, kExtHi, C_ext, C_ext_err, chi2ndf_ext), "l");
    leg->AddEntry(l_fit_uw,
        Form("pol0 UW  [%.2f, %.2f] unweighted: C_{UW} = %.3f #pm %.3f, "
             "#chi^{2}/ndf = %.1f/%d",
             kExtLo, kExtHi, C_uw, C_uw_err, chi2_uw, ndf_uw), "l");
    leg->AddEntry((TObject*)nullptr,
        Form("#DeltaC = C_{EXT} - C_{MID} = %+.3f  (%.2f #sigma)",
             dC, dC_signif), "");
    leg->Draw();

    // Annotation: per-region means.
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.026);
    tx.SetTextColor(kBlack);
    tx.DrawLatex(0.14, 0.84,
        Form("LOW  [0.00, %.2f]: <r> = %.3f #pm %.3f",
             kLowHi, mean_low, sig_low));
    tx.DrawLatex(0.14, 0.80,
        Form("MID  [%.2f, %.2f]: <r> = %.3f #pm %.3f",
             kLowHi, kHighLo, mean_mid, sig_mid));
    tx.DrawLatex(0.14, 0.76,
        Form("HIGH [%.2f, %.2f]: <r> = %.3f #pm %.3f",
             kHighLo, kXmax, mean_high, sig_high));

    // ---------- Save ----------
    const std::string base = "plots/output/mass_ee_pt_ratio_innermdc";
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
    std::cout << "\nSaved: " << base << ".{pdf,png}\n";

    f_em->Close();
    f_pp->Close();
    f_mm->Close();
}
