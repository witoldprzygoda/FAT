// mass_ee_signal_trigger_corrected_step5pts.C
// =========================================================================
// EXP analysis: Fine-binned (280 x 5 MeV) PT3 CB-subtracted m_ee signal,
// drawn alongside the SAME spectrum with a trigger-bias correction applied
// bin-by-bin. Mirror of mass_ee_signal_trigger_corrected_step5.C, BUT the
// correction comes from MEASURED step5 SIM PT2/PT3 RATIO POINTS (cascade
// v5 220-bin variable layout) scaled by a SINGLE factor k_simpts to match
// the EXP signal ratio integral over [0.15, 0.70].
//
//   For each fine PT3 bin (280 x 5 MeV on [0, 1.4]) with center x:
//       correction(x) = k_simpts * sim_ratio_at(x)
//          (nearest-neighbour lookup from sim ratio TH1, NO LOW/HIGH split)
//       corrected_signal(x) = signal(x) * correction(x)
//       sigma_corrected(x)  = correction(x) * sigma_signal(x)
//
//   k_simpts = I_expRatio[0.15, 0.70] / I_simpoints[0.15, 0.70]
//
// where I_expRatio is the 58-bin EXP signal ratio integral and
// I_simpoints is the 220-bin sim ratio integral, both via
// Σ (bin_content · bin_width) over bins fully contained in [0.15, 0.70].
//
// Output:
//   plots/output/mass_ee_signal_trigger_corrected_step5pts.{pdf,png}
//
// Usage:
//   root -l -b -q plots/mass_ee_signal_trigger_corrected_step5pts.C
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

// SIM inputs.
constexpr const char* kFileSim    = "output_epem_sim.root";
constexpr const char* kTreeSimRec = "dilepton_nt";

// Step5 (same_vertex) cumulative purity gate.
constexpr const char* kSimStep5Cut =
    "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)*"
    "(epem_same_vertex==1)";

constexpr const char* kCutPT3 = "trigbit==8192";
constexpr const char* kCutPT2 = "trigbit==4096";

constexpr double kTrigCorr = 63.0;

// -------- Fine binning: 280 uniform bins x 5 MeV over [0, 1.4] --------
constexpr int    kFineNb   = 280;
constexpr double kFineXmin = 0.0;
constexpr double kFineXmax = 1.4;

// -------- Sparse binning: 58 variable bins (matches signal ratio macro) -
constexpr int    kSparseNb   = 58;
constexpr double kSparseXmin = 0.0;
constexpr double kSparseXmax = 1.4;

// SIM variable binning (cascade v5: 220 bins on [0, 1.4]).
constexpr int kSimNb = 220;
const double* getSimEdges() {
    static double edges[kSimNb + 1];
    static bool   built = false;
    if (!built) {
        for (int i = 0; i <= 60;  ++i) edges[i]       = i * 0.0025;
        for (int i = 1; i <= 140; ++i) edges[60 + i]  = 0.15 + i * 0.005;
        for (int i = 1; i <= 12;  ++i) edges[200 + i] = 0.85 + i * (0.25 / 12.0);
        for (int i = 1; i <= 8;   ++i) edges[212 + i] = 1.10 + i * (0.30 / 8.0);
        built = true;
    }
    return edges;
}
constexpr double kSimMinDen = 5.0;

// LOW/HIGH split for integral reporting only (correction has NO split).
constexpr double kSplitX = 0.15;

// Integral-matching range for k_simpts.
constexpr double kIntLo = 0.15, kIntHi = 0.70;

// Build the 58-bin variable-width edge array (matches signal ratio macro).
const double* getSparseEdges() {
    static double edges[kSparseNb + 1];
    static bool   built = false;
    if (!built) {
        int k = 0;
        for (int i = 0; i <= 20; ++i) edges[k++] = i * 0.0025;
        const double w = (1.4 - 0.05) / 38.0;
        for (int i = 1; i <= 38; ++i) edges[k++] = 0.05 + i * w;
        built = true;
    }
    return edges;
}

// Draw m_ee with the given cut into a histogram of the specified layout.
TH1D* drawFine(TTree* t, const char* cut, const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kFineNb, kFineXmin, kFineXmax);
    h->Sumw2();
    t->Draw(("m_ee>>" + name).c_str(), cut, "goff");
    h->SetDirectory(nullptr);
    return h;
}
TH1D* drawSparse(TTree* t, const char* cut, const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kSparseNb, getSparseEdges());
    h->Sumw2();
    t->Draw(("m_ee>>" + name).c_str(), cut, "goff");
    h->SetDirectory(nullptr);
    return h;
}

// CB = 2*sqrt(N++ * N--) per bin with proper error propagation.
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

// Signal = epem - CB with quadrature error propagation.
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

// Bin-by-bin ratio = factor * num / den with relative-error propagation.
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
        if (d == 0.0 || !std::isfinite(d) || !std::isfinite(n)) continue;
        const double val = factor * n / d;
        const double rel = std::sqrt(
            (n != 0.0 ? std::pow(en / n, 2) : 0.0) +
            std::pow(ed / d, 2));
        r->SetBinContent(b, val);
        r->SetBinError  (b, std::abs(val) * rel);
    }
    return r;
}

// Divide each bin (content & error) by its bin width.
void normalizeByBinWidth(TH1D* h) {
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double w = h->GetBinWidth(b);
        if (w > 0.0) {
            h->SetBinContent(b, h->GetBinContent(b) / w);
            h->SetBinError  (b, h->GetBinError(b)   / w);
        }
    }
}

void styleHist(TH1D* h, Color_t c, Style_t m) {
    h->SetLineColor(c);
    h->SetMarkerColor(c);
    h->SetMarkerStyle(m);
    h->SetMarkerSize(0.7);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.040);
    h->GetYaxis()->SetLabelSize(0.040);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

// Sum bin contents over [xlo, xhi] (using bin centers to decide inclusion).
double sumByCenter(TH1D* h, double xlo, double xhi) {
    double s = 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double xc = h->GetXaxis()->GetBinCenter(b);
        if (xc >= xlo && xc <= xhi) s += h->GetBinContent(b);
    }
    return s;
}

// Integral helper: sum (bin_content * bin_width) over bins fully contained
// in [xlo, xhi].
double histIntegralXBinWidth(TH1D* h, double xlo, double xhi) {
    double s = 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double lo = h->GetXaxis()->GetBinLowEdge(b);
        const double hi = h->GetXaxis()->GetBinUpEdge (b);
        if (lo < xlo - 1e-9) continue;
        if (hi > xhi + 1e-9) continue;
        s += h->GetBinContent(b) * h->GetBinWidth(b);
    }
    return s;
}

// Build sim ratio histogram on cascade-v5 220-bin layout.
TH1D* buildSimFullRatio() {
    TFile* f_sim = TFile::Open(kFileSim, "READ");
    if (!f_sim || f_sim->IsZombie()) {
        std::cerr << "[sim] Cannot open sim file: " << kFileSim << "\n";
        if (f_sim) f_sim->Close();
        return nullptr;
    }
    TTree* t_sim = dynamic_cast<TTree*>(f_sim->Get(kTreeSimRec));
    if (!t_sim) {
        std::cerr << "[sim] Tree '" << kTreeSimRec
                  << "' missing in sim file\n";
        f_sim->Close();
        return nullptr;
    }
    TH1D* h_sim_PT3 = new TH1D("h_sim_full_PT3_corr", "",
                               kSimNb, getSimEdges());
    TH1D* h_sim_PT2 = new TH1D("h_sim_full_PT2_corr", "",
                               kSimNb, getSimEdges());
    h_sim_PT3->Sumw2();
    h_sim_PT2->Sumw2();

    const std::string w_PT3 =
        std::string("(pt3==1)*") + kSimStep5Cut + "*sim_genweight";
    const std::string w_PT2 =
        std::string("(pt2==1)*") + kSimStep5Cut + "*sim_genweight";
    t_sim->Draw("m_ee>>h_sim_full_PT3_corr", w_PT3.c_str(), "goff");
    t_sim->Draw("m_ee>>h_sim_full_PT2_corr", w_PT2.c_str(), "goff");
    h_sim_PT3->SetDirectory(nullptr);
    h_sim_PT2->SetDirectory(nullptr);

    TH1D* r = static_cast<TH1D*>(h_sim_PT3->Clone("h_sim_full_ratio_corr"));
    r->SetDirectory(nullptr);
    r->Reset();
    int n_drawn = 0;
    for (int b = 1; b <= kSimNb; ++b) {
        const double n2 = h_sim_PT2->GetBinContent(b);
        const double e2 = h_sim_PT2->GetBinError  (b);
        const double n3 = h_sim_PT3->GetBinContent(b);
        const double e3 = h_sim_PT3->GetBinError  (b);
        if (n3 < kSimMinDen || !std::isfinite(n3) || !std::isfinite(n2)) {
            r->SetBinContent(b, 0.0);
            r->SetBinError  (b, 0.0);
            continue;
        }
        const double val = n2 / n3;
        const double rel = std::sqrt(
            (n2 != 0.0 ? std::pow(e2 / n2, 2) : 0.0) +
            std::pow(e3 / n3, 2));
        r->SetBinContent(b, val);
        r->SetBinError  (b, std::abs(val) * rel);
        ++n_drawn;
    }
    printf("[sim] FULL 220 bins drawn (PT3>=%.0f): %d / %d\n",
           kSimMinDen, n_drawn, kSimNb);
    delete h_sim_PT3;
    delete h_sim_PT2;
    f_sim->Close();
    return r;
}

// Nearest-neighbour lookup with outward walk for missing bins.
double simRatioAt(TH1D* h_sim_ratio, double x) {
    if (!h_sim_ratio) return 0.0;
    int b = h_sim_ratio->FindBin(x);
    if (b < 1) b = 1;
    if (b > h_sim_ratio->GetNbinsX()) b = h_sim_ratio->GetNbinsX();
    if (h_sim_ratio->GetBinContent(b) > 0.0) {
        return h_sim_ratio->GetBinContent(b);
    }
    const int nb = h_sim_ratio->GetNbinsX();
    for (int d = 1; d <= nb; ++d) {
        const int bl = b - d;
        const int bu = b + d;
        if (bl >= 1 && h_sim_ratio->GetBinContent(bl) > 0.0) {
            return h_sim_ratio->GetBinContent(bl);
        }
        if (bu <= nb && h_sim_ratio->GetBinContent(bu) > 0.0) {
            return h_sim_ratio->GetBinContent(bu);
        }
    }
    return 0.0;
}

}  // anonymous namespace

void mass_ee_signal_trigger_corrected_step5pts() {
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
    std::cout << "  epem : " << kFileEm
              << " (" << t_em->GetEntries() << " entries)\n";
    std::cout << "  epep : " << kFilePp
              << " (" << t_pp->GetEntries() << " entries)\n";
    std::cout << "  emem : " << kFileMm
              << " (" << t_mm->GetEntries() << " entries)\n";

    std::cout << "Cuts: PT3='" << kCutPT3
              << "'  PT2='" << kCutPT2 << "'\n";
    std::cout << "SIM step5 (same_vertex) cut + sim_genweight applied "
                 "to sim ratio.\n";
    std::cout << "Correction: corr(x) = k_simpts * sim_ratio_at(x)  "
                 "[nearest-neighbour lookup, full mass range, NO LOW/HIGH "
                 "split].\n\n";

    gSystem->mkdir("plots/output", true);

    // =====================================================================
    // SECTION A: Fine-binned PT3 signal mass spectrum (280 x 5 MeV)
    // =====================================================================
    std::cout << "============ SECTION A: fine-binned PT3 signal ============\n";

    TH1D* h_pt3_epem = drawFine(t_em, kCutPT3, "h_pt3_epem");
    TH1D* h_pt3_epep = drawFine(t_pp, kCutPT3, "h_pt3_epep");
    TH1D* h_pt3_emem = drawFine(t_mm, kCutPT3, "h_pt3_emem");

    TH1D* h_pt3_cb     = makeCB(h_pt3_epep, h_pt3_emem, "h_pt3_cb");
    TH1D* h_pt3_signal = makeSignal(h_pt3_epem, h_pt3_cb, "h_pt3_signal");

    std::cout << "Fine PT3 raw yields:\n";
    std::cout << "  epem = " << h_pt3_epem->Integral() << "\n";
    std::cout << "  epep = " << h_pt3_epep->Integral() << "\n";
    std::cout << "  emem = " << h_pt3_emem->Integral() << "\n";
    std::cout << "  CB     = " << h_pt3_cb->Integral() << "\n";
    std::cout << "  signal = " << h_pt3_signal->Integral() << "\n";

    // =====================================================================
    // SECTION B: Build EXP signal ratio (58-bin) and SIM ratio (220-bin)
    //            compute k_simpts via integral matching
    // =====================================================================
    std::cout << "\n========== SECTION B: build ratios, compute k_simpts ==========\n";

    TH1D* h_sp_em_PT3 = drawSparse(t_em, kCutPT3, "h_sp_em_PT3");
    TH1D* h_sp_em_PT2 = drawSparse(t_em, kCutPT2, "h_sp_em_PT2");
    TH1D* h_sp_pp_PT3 = drawSparse(t_pp, kCutPT3, "h_sp_pp_PT3");
    TH1D* h_sp_pp_PT2 = drawSparse(t_pp, kCutPT2, "h_sp_pp_PT2");
    TH1D* h_sp_mm_PT3 = drawSparse(t_mm, kCutPT3, "h_sp_mm_PT3");
    TH1D* h_sp_mm_PT2 = drawSparse(t_mm, kCutPT2, "h_sp_mm_PT2");

    TH1D* h_sp_cb_PT3  = makeCB(h_sp_pp_PT3, h_sp_mm_PT3, "h_sp_cb_PT3");
    TH1D* h_sp_cb_PT2  = makeCB(h_sp_pp_PT2, h_sp_mm_PT2, "h_sp_cb_PT2");
    TH1D* h_sp_sig_PT3 = makeSignal(h_sp_em_PT3, h_sp_cb_PT3, "h_sp_sig_PT3");
    TH1D* h_sp_sig_PT2 = makeSignal(h_sp_em_PT2, h_sp_cb_PT2, "h_sp_sig_PT2");

    // ratio(b) = 63 * sig_PT2(b) / sig_PT3(b)  (signal ratio, EXP)
    TH1D* h_sp_ratio = makeRatio(h_sp_sig_PT2, h_sp_sig_PT3, kTrigCorr,
                                 "h_sp_ratio");

    // Build sim ratio (220-bin).
    TH1D* h_sim_ratio = buildSimFullRatio();
    if (!h_sim_ratio) {
        std::cerr << "Sim ratio could not be built; aborting.\n";
        return;
    }

    // ---------- INTEGRAL-MATCHED k_simpts over [kIntLo, kIntHi] ----------
    const double I_expRatio_mid =
        histIntegralXBinWidth(h_sp_ratio, kIntLo, kIntHi);
    const double I_simpoints_mid =
        histIntegralXBinWidth(h_sim_ratio, kIntLo, kIntHi);
    const double k_simpts = (I_simpoints_mid != 0.0)
        ? I_expRatio_mid / I_simpoints_mid
        : 0.0;

    printf("\nINTEGRAL-MATCHED k_simpts (range [%.2f, %.2f]):\n",
           kIntLo, kIntHi);
    printf("  I_expRatio  (sum exp_signal_ratio(b) x binwidth) = %.6f\n",
           I_expRatio_mid);
    printf("  I_simpoints (sum sim_ratio(b) x binwidth)        = %.6f\n",
           I_simpoints_mid);
    printf("  k_simpts    = I_expRatio / I_simpoints            = %.6f\n",
           k_simpts);

    // =====================================================================
    // SECTION C: Apply POINTS-based correction to fine-binned PT3 signal
    //            (no LOW/HIGH split — single sim-points formula everywhere)
    // =====================================================================
    std::cout << "\n========== SECTION C: apply POINTS correction ==========\n";

    TH1D* h_pt3_signal_corrected = static_cast<TH1D*>(
        h_pt3_signal->Clone("h_pt3_signal_corrected"));
    h_pt3_signal_corrected->SetDirectory(nullptr);
    h_pt3_signal_corrected->Reset();

    int n_filled = 0;
    int n_missing = 0;
    for (int b = 1; b <= h_pt3_signal->GetNbinsX(); ++b) {
        const double xc  = h_pt3_signal->GetXaxis()->GetBinCenter(b);
        const double sv  = h_pt3_signal->GetBinContent(b);
        const double sev = h_pt3_signal->GetBinError  (b);
        const double sr  = simRatioAt(h_sim_ratio, xc);
        const double corr = k_simpts * sr;
        if (corr <= 0.0 || !std::isfinite(corr)) {
            h_pt3_signal_corrected->SetBinContent(b, 0.0);
            h_pt3_signal_corrected->SetBinError  (b, 0.0);
            ++n_missing;
            continue;
        }
        const double val = sv * corr;
        const double err = std::abs(corr) * sev;   // sim error neglected
        h_pt3_signal_corrected->SetBinContent(b, val);
        h_pt3_signal_corrected->SetBinError  (b, err);
        ++n_filled;
    }
    printf("Correction applied: %d bins filled, %d bins missing (corr=0)\n",
           n_filled, n_missing);

    // =====================================================================
    // SECTION E: Integrals
    // =====================================================================
    std::cout << "\n========== SECTION E: integrals ==========\n";
    const double I_raw       = h_pt3_signal->Integral();
    const double I_corrected = h_pt3_signal_corrected->Integral();
    const double delta       = I_corrected - I_raw;
    const double ratio       = (I_raw != 0.0) ? I_corrected / I_raw : 0.0;
    printf("I_raw       = %.6f\n", I_raw);
    printf("I_corrected = %.6f\n", I_corrected);
    printf("delta       = I_corrected - I_raw = %.6f\n", delta);
    if (ratio > 1.0) {
        printf("ratio       = I_corrected / I_raw = %.6f  -> corrected is "
               "%.3f%% higher\n",
               ratio, 100.0 * (ratio - 1.0));
    } else if (ratio < 1.0 && ratio > 0.0) {
        printf("ratio       = I_corrected / I_raw = %.6f  -> corrected is "
               "%.3f%% lower\n",
               ratio, 100.0 * (1.0 - ratio));
    } else {
        printf("ratio       = I_corrected / I_raw = %.6f\n", ratio);
    }

    const double xlo = kFineXmin + 1e-9;
    const double xhi = kFineXmax - 1e-9;
    const double I_raw_LOW  = sumByCenter(h_pt3_signal, xlo, kSplitX);
    const double I_raw_HIGH = sumByCenter(h_pt3_signal, kSplitX, xhi);
    const double I_cor_LOW  = sumByCenter(h_pt3_signal_corrected, xlo, kSplitX);
    const double I_cor_HIGH = sumByCenter(h_pt3_signal_corrected, kSplitX, xhi);
    printf("\nSplit at %.2f:\n", kSplitX);
    printf("  I_raw_LOW  = %.6f   I_cor_LOW  = %.6f   (delta %+.6f, ratio %.6f)\n",
           I_raw_LOW, I_cor_LOW, I_cor_LOW - I_raw_LOW,
           (I_raw_LOW != 0.0 ? I_cor_LOW / I_raw_LOW : 0.0));
    printf("  I_raw_HIGH = %.6f   I_cor_HIGH = %.6f   (delta %+.6f, ratio %.6f)\n",
           I_raw_HIGH, I_cor_HIGH, I_cor_HIGH - I_raw_HIGH,
           (I_raw_HIGH != 0.0 ? I_cor_HIGH / I_raw_HIGH : 0.0));

    // =====================================================================
    // SECTION D: Drawing
    // =====================================================================
    std::cout << "\n========== SECTION D: drawing ==========\n";

    TH1D* h_pt3_signal_w_norm = static_cast<TH1D*>(
        h_pt3_signal->Clone("h_pt3_signal_w_norm"));
    h_pt3_signal_w_norm->SetDirectory(nullptr);
    TH1D* h_pt3_signal_corrected_w_norm = static_cast<TH1D*>(
        h_pt3_signal_corrected->Clone("h_pt3_signal_corrected_w_norm"));
    h_pt3_signal_corrected_w_norm->SetDirectory(nullptr);
    normalizeByBinWidth(h_pt3_signal_w_norm);
    normalizeByBinWidth(h_pt3_signal_corrected_w_norm);

    styleHist(h_pt3_signal_w_norm,           kBlack, 20);
    styleHist(h_pt3_signal_corrected_w_norm, kRed,   20);

    TCanvas* c = new TCanvas("c_trigcorr_step5pts",
                             "PT3 signal raw vs trigger-bias corrected (step5 POINTS)",
                             1100, 700);
    c->SetMargin(0.13, 0.04, 0.13, 0.10);
    gPad->SetLogy(true);

    double y_max_raw = h_pt3_signal_w_norm->GetMaximum();
    double y_max_cor = h_pt3_signal_corrected_w_norm->GetMaximum();
    double y_max     = 1.5 * std::max(y_max_raw, y_max_cor);
    double y_min_pos = 1e300;
    for (int b = 1; b <= h_pt3_signal_w_norm->GetNbinsX(); ++b) {
        const double v1 = h_pt3_signal_w_norm->GetBinContent(b);
        const double v2 = h_pt3_signal_corrected_w_norm->GetBinContent(b);
        if (v1 > 0.0 && v1 < y_min_pos) y_min_pos = v1;
        if (v2 > 0.0 && v2 < y_min_pos) y_min_pos = v2;
    }
    if (!std::isfinite(y_min_pos) || y_min_pos > 1e299) y_min_pos = 0.1;
    double y_min = std::max(0.1, 0.3 * y_min_pos);

    h_pt3_signal_w_norm->SetTitle(
        "m_{ee} PT3 signal - raw vs trigger-bias corrected (step5 POINTS);"
        "m_{ee} [GeV/c^{2}];weighted entries / GeV/c^{2}");
    h_pt3_signal_w_norm->GetXaxis()->SetRangeUser(kFineXmin, kFineXmax);
    h_pt3_signal_w_norm->GetYaxis()->SetRangeUser(y_min, y_max);
    h_pt3_signal_w_norm->Draw("E1");
    h_pt3_signal_corrected_w_norm->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.45, 0.72, 0.95, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.028);
    leg->AddEntry(h_pt3_signal_w_norm,
                  "PT3 signal (raw)", "lpe");
    leg->AddEntry(h_pt3_signal_corrected_w_norm,
                  Form("PT3 trigger bias corrected "
                       "(sim step5 POINTS, k_{simpts} = %.3f)", k_simpts),
                  "lpe");
    leg->Draw();

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.026);
    tx.DrawLatex(0.16, 0.86,
        Form("k_{simpts} (norm: integral [%.2f, %.2f]) = %.6f",
             kIntLo, kIntHi, k_simpts));
    tx.DrawLatex(0.16, 0.82,
        Form("I_{exp}^{ratio} = %.4f,  I_{sim}^{points} = %.4f",
             I_expRatio_mid, I_simpoints_mid));
    tx.DrawLatex(0.16, 0.78,
        Form("I_{raw} = %.0f,  I_{cor} = %.0f,  ratio = %.4f",
             I_raw, I_corrected, ratio));
    tx.DrawLatex(0.16, 0.74,
        Form("LOW [0, %.2f]: I_{raw}=%.0f, I_{cor}=%.0f",
             kSplitX, I_raw_LOW, I_cor_LOW));
    tx.DrawLatex(0.16, 0.70,
        Form("HIGH (%.2f, 1.4]: I_{raw}=%.0f, I_{cor}=%.0f",
             kSplitX, I_raw_HIGH, I_cor_HIGH));

    // ---------- Save ----------
    const std::string out_base =
        "plots/output/mass_ee_signal_trigger_corrected_step5pts";
    c->SaveAs((out_base + ".pdf").c_str());
    c->SaveAs((out_base + ".png").c_str());
    std::cout << "Saved: " << out_base << ".{pdf,png}\n";

    f_em->Close();
    f_pp->Close();
    f_mm->Close();
}
