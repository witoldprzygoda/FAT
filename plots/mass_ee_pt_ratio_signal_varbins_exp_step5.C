// mass_ee_pt_ratio_signal_varbins_exp_step5.C
// =========================================================================
// EXP analysis: 63 · N_PT2 / N_PT3 ratio with variable binning. Mirror of
// mass_ee_pt_ratio_signal_varbins_exp.C using step5 (same_vertex) sim fit
// constants and INTEGRAL-MATCHING normalization over [0.2, 0.8] instead of
// chi² fit over [0.10, 1.40].
//
// Output:
//   plots/output/mass_ee_pt_ratio_varbins_exp_signal_step5.{pdf,png}
//
// Step5 (same_vertex) sim cut for the LOW red-circle overlay:
//   "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)*"
//   "(epem_same_vertex==1)"
//
// Usage: root -l -b -q plots/mass_ee_pt_ratio_signal_varbins_exp_step5.C
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
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr const char* kFileEm = "output_epem_exp.root";
constexpr const char* kFilePp = "output_epep_exp.root";
constexpr const char* kFileMm = "output_emem_exp.root";
constexpr const char* kTree   = "dilepton_nt";

// SIM overlay inputs — used to draw a subtle reference in the LOW [0.0, 0.15]
// region of the right (ratio) panel.
constexpr const char* kFileSim    = "output_epem_sim.root";
constexpr const char* kTreeSimRec = "dilepton_nt";

// Step5 (same_vertex) cumulative purity gate.
constexpr const char* kSimStep5Cut =
    "(eVertReco_z>-500)*(isBest==1)*(ep_sim_id==2)*(em_sim_id==3)*"
    "(epem_same_vertex==1)";

// SIM low-region binning: 30 bins × 5 MeV over [0.00, 0.15].
constexpr int    kSimLowNb   = 30;
constexpr double kSimLowXmin = 0.00;
constexpr double kSimLowXmax = 0.15;
// Minimum sim PT3 counts per bin to draw an overlay point.
constexpr double kSimMinDen  = 5.0;

constexpr int    kNb   = 58;
constexpr double kXmin = 0.0;
constexpr double kXmax = 1.4;

// LOW-normalized fit range — [0.00, 0.10].
constexpr double kLowLo = 0.00, kLowHi = 0.10;

// Trigger downscale: PT2 is downscaled by 64, so we multiply by 63.
constexpr double kTrigCorr = 63.0;

// Trigger bit cuts (EXP) — NO additional selection.
constexpr const char* kCutPT3 = "trigbit==8192";
constexpr const char* kCutPT2 = "trigbit==4096";

// Native MID fit range — [0.10, 0.80].
constexpr double kMdLo = 0.10, kMdHi = 0.80;

// Integral-matching range for k_step5_mid (NEW).
// Match k via integral over the MID region (where sim MID was originally
// fitted in plots/mass_ee_pt_ratio_purity_cascade_sim_3piece_v5.C).
// The SAME k is then applied to HIGH part (HIGH becomes a sim-predicted
// extrapolation, not separately re-normalized).
constexpr double kIntLo = 0.15, kIntHi = 0.70;

constexpr Color_t kColPT3   = kBlack;
constexpr Color_t kColPT2   = kRed   + 1;
constexpr Color_t kColRatio = kBlue  + 1;

// Native MID free quad colour & style.
const Color_t kColMid       = kAzure  + 1;   // DASHED
// Sim recipe display piece colours (two separate TF1s).
const Color_t kColSimMid    = kAzure  + 3;   // SOLID, MID portion of sim shape
const Color_t kColSimHigh   = kOrange + 8;   // SOLID, HIGH portion of sim shape

// -------- SCALED SIM RECIPE constants (from sim cascade v5 step5 REC) --------
// step5 = same_vertex cut applied on top of step4.
constexpr double SIM_JOIN_X  = 0.70;
constexpr double SIM_HI_A    =  1.467;
constexpr double SIM_HI_B    = -0.803;
constexpr double SIM_HI_C    =  0.738;
constexpr double SIM_MID_B   =  0.806;
constexpr double SIM_MID_C   =  1.444;
// sim_join_y (HIGH evaluated at SIM_JOIN_X):
const double SIM_JOIN_Y =
    SIM_HI_A + SIM_HI_B * SIM_JOIN_X + SIM_HI_C * SIM_JOIN_X * SIM_JOIN_X;

// Display-only TF1 functions (one parameter k each, evaluated on a fixed
// sim-piece formula).
double simMidDisplay(double* x, double* par) {
    const double t = x[0];
    const double k = par[0];
    const double d = t - SIM_JOIN_X;
    return k * (SIM_JOIN_Y + SIM_MID_B * d + SIM_MID_C * d * d);
}
double simHighDisplay(double* x, double* par) {
    const double t = x[0];
    const double k = par[0];
    return k * (SIM_HI_A + SIM_HI_B * t + SIM_HI_C * t * t);
}

// Bare piecewise sim shape (no scale factor) — used by the LOW-normalized
// refit and to predict 63*PT2_hypo from PT3 on the left panel.
double simShape(double x) {
    if (x < SIM_JOIN_X) {
        const double d = x - SIM_JOIN_X;
        return SIM_JOIN_Y + SIM_MID_B * d + SIM_MID_C * d * d;
    }
    return SIM_HI_A + SIM_HI_B * x + SIM_HI_C * x * x;
}

// TF1-compatible: y = k * simShape(x), where k = par[0].
double simShapeScaled(double* x, double* par) {
    return par[0] * simShape(x[0]);
}

// Closed-form integral of simShape(x) over [a, b]. Piecewise:
//   for x < SIM_JOIN_X : f(x) = SIM_JOIN_Y + SIM_MID_B*(x-J) + SIM_MID_C*(x-J)^2
//   for x >= SIM_JOIN_X: f(x) = SIM_HI_A   + SIM_HI_B*x      + SIM_HI_C*x^2
// Antiderivative on MID with substitution u = x - J:
//   F_mid(x) = SIM_JOIN_Y*x + 0.5*SIM_MID_B*(x-J)^2 + (1/3)*SIM_MID_C*(x-J)^3
// Antiderivative on HIGH:
//   F_hi(x)  = SIM_HI_A*x + 0.5*SIM_HI_B*x^2 + (1/3)*SIM_HI_C*x^3
double simShapeAntiMid(double x) {
    const double d = x - SIM_JOIN_X;
    return SIM_JOIN_Y * x
         + 0.5 * SIM_MID_B * d * d
         + (1.0 / 3.0) * SIM_MID_C * d * d * d;
}
double simShapeAntiHigh(double x) {
    return SIM_HI_A * x
         + 0.5 * SIM_HI_B * x * x
         + (1.0 / 3.0) * SIM_HI_C * x * x * x;
}
double simShapeIntegral(double a, double b) {
    if (b <= a) return 0.0;
    double I = 0.0;
    const double J = SIM_JOIN_X;
    if (b <= J) {
        I = simShapeAntiMid(b) - simShapeAntiMid(a);
    } else if (a >= J) {
        I = simShapeAntiHigh(b) - simShapeAntiHigh(a);
    } else {
        I  = simShapeAntiMid(J)  - simShapeAntiMid(a);
        I += simShapeAntiHigh(b) - simShapeAntiHigh(J);
    }
    return I;
}

// Build the variable-width bin edge array (58 bins total).
//   [0.00, 0.05] : 20 bins of 2.5 MeV       (DENSE LOW — restricted)
//   [0.05, 1.40] : 38 bins of ~35.526 MeV   (UNIFORM, MID/HIGH-like width)
const double* getEdges() {
    static double edges[kNb + 1];
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

// Fill a histogram from a TTree using TTree::Draw with cut. Sumw2 first.
TH1D* drawWithCut(TTree* t, const char* expr,
                  const char* cut, const std::string& name) {
    TH1D* h = new TH1D(name.c_str(), "", kNb, getEdges());
    h->Sumw2();
    t->Draw((std::string(expr) + ">>" + name).c_str(), cut, "goff");
    h->SetDirectory(nullptr);
    return h;
}

// Compute CB_t = 2*sqrt(N_++ · N_--) per bin with proper error propagation.
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

// Signal = epem − CB with quadrature errors.
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

// Bin-by-bin ratio = factor · num / den with relative-error propagation.
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

// Divide each bin (content & error) by its bin width — spectra display.
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
    h->SetMarkerSize(0.9);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

// Free quadratic [0]+[1]*x+[2]*x*x on [xlo, xhi].
struct QuadResult {
    double a, a_err;
    double b, b_err;
    double c, c_err;
    double chi2_ndf;
    int    ndf;
    TF1*   f;
};

QuadResult fitQuadFree(TH1D* h, double xlo, double xhi,
                       Color_t color, Style_t style, const char* tag) {
    const std::string fname = std::string(h->GetName()) + "_fq_" + tag;
    TF1* f = new TF1(fname.c_str(), "[0]+[1]*x+[2]*x*x", xlo, xhi);
    f->SetLineColor(color);
    f->SetLineWidth(3);
    f->SetLineStyle(style);
    h->Fit(f, "RQ0");
    QuadResult r;
    r.a        = f->GetParameter(0);
    r.a_err    = f->GetParError(0);
    r.b        = f->GetParameter(1);
    r.b_err    = f->GetParError(1);
    r.c        = f->GetParameter(2);
    r.c_err    = f->GetParError(2);
    r.ndf      = f->GetNDF();
    r.chi2_ndf = (r.ndf > 0) ? f->GetChisquare() / r.ndf : 0.0;
    r.f        = f;
    return r;
}

// Sum of bin content within [xlo, xhi] (inclusive on both ends).
double sumInRange(TH1D* h, double xlo, double xhi) {
    double s = 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double lo = h->GetXaxis()->GetBinLowEdge(b);
        const double hi = h->GetXaxis()->GetBinUpEdge(b);
        if (lo >= xlo - 1e-9 && hi <= xhi + 1e-9) {
            s += h->GetBinContent(b);
        }
    }
    return s;
}

// Integral helper: sum over hist bins fully contained in [xlo, xhi] of
// (bin content * bin width). Useful for normalising a step-function "shape"
// embodied by a histogram.
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

// Print per-bin diagnostics in [xlo, xhi]: for every bin whose edges fall
// within the range, dump (low_edge, high_edge, content).
void dumpBinsInRange(const char* label, TH1D* h, double xlo, double xhi,
                     double mult = 1.0) {
    printf("  %s per-bin in [%.2f, %.2f]:\n", label, xlo, xhi);
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double lo = h->GetXaxis()->GetBinLowEdge(b);
        const double hi = h->GetXaxis()->GetBinUpEdge(b);
        if (lo >= xlo - 1e-9 && hi <= xhi + 1e-9) {
            const double c = h->GetBinContent(b) * mult;
            printf("    [%.4f, %.4f] -> %.3f\n", lo, hi, c);
        }
    }
}

// ---------------------------------------------------------------------------
// Build a SIM low-region ratio histogram for overlay on the right panel.
// Uses step5 (same_vertex) sim cut.
// ---------------------------------------------------------------------------
TH1D* buildSimLowRatio(const std::string& tag) {
    TFile* f_sim = TFile::Open(kFileSim, "READ");
    if (!f_sim || f_sim->IsZombie()) {
        std::cerr << "[sim overlay] Cannot open sim file: "
                  << kFileSim << "\n";
        if (f_sim) f_sim->Close();
        return nullptr;
    }
    TTree* t_sim = dynamic_cast<TTree*>(f_sim->Get(kTreeSimRec));
    if (!t_sim) {
        std::cerr << "[sim overlay] Tree '" << kTreeSimRec
                  << "' missing in sim file\n";
        f_sim->Close();
        return nullptr;
    }

    // Build uniform 30-bin edges over [0.00, 0.15] (5 MeV bins).
    double edges[kSimLowNb + 1];
    for (int i = 0; i <= kSimLowNb; ++i) {
        edges[i] = kSimLowXmin +
                   (kSimLowXmax - kSimLowXmin) * static_cast<double>(i) /
                   static_cast<double>(kSimLowNb);
    }

    const std::string nm3 = "h_sim_low_PT3_" + tag;
    const std::string nm2 = "h_sim_low_PT2_" + tag;
    TH1D* h_sim_PT3 = new TH1D(nm3.c_str(), "", kSimLowNb, edges);
    TH1D* h_sim_PT2 = new TH1D(nm2.c_str(), "", kSimLowNb, edges);
    h_sim_PT3->Sumw2();
    h_sim_PT2->Sumw2();

    // Apply step5 cumulative cut + trigger gate + sim_genweight, project m_ee.
    const std::string w_PT3 =
        std::string("(pt3==1)*") + kSimStep5Cut + "*sim_genweight";
    const std::string w_PT2 =
        std::string("(pt2==1)*") + kSimStep5Cut + "*sim_genweight";
    t_sim->Draw(("m_ee>>" + nm3).c_str(), w_PT3.c_str(), "goff");
    t_sim->Draw(("m_ee>>" + nm2).c_str(), w_PT2.c_str(), "goff");
    h_sim_PT3->SetDirectory(nullptr);
    h_sim_PT2->SetDirectory(nullptr);

    // Compute ratio = PT2 / PT3 with independent-counts error propagation.
    TH1D* r = static_cast<TH1D*>(
        h_sim_PT3->Clone(("h_sim_low_ratio_" + tag).c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    int n_drawn = 0;
    for (int b = 1; b <= kSimLowNb; ++b) {
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

    printf("[sim overlay] LOW [0.00, 0.15] 30x5 MeV bins drawn: %d / %d\n",
           n_drawn, kSimLowNb);

    delete h_sim_PT3;
    delete h_sim_PT2;
    f_sim->Close();
    return r;
}

// ---------------------------------------------------------------------------
// Produce one full output file (PDF + PNG) for the given numerator/denominator
// histograms and the corresponding display spectra.
// ---------------------------------------------------------------------------
void makeOnePlot(TH1D* h_num, TH1D* h_den,
                 TH1D* h_disp_PT3, TH1D* h_disp_PT2,
                 const char* left_title,
                 const char* right_title,
                 const char* right_ylabel,
                 const std::string& base_path,
                 const std::string& tag,
                 const std::string& variant_label) {

    std::cout << "\n========== " << variant_label << " ==========\n";

    // ---------- Ratio = 63 · num / den (before width-norm) ----------
    TH1D* r_ratio = makeRatio(h_num, h_den, kTrigCorr,
                              "r_ratio_" + tag);
    styleHist(r_ratio, kColRatio, 20);

    // ---------- Native MID free quad ----------
    QuadResult mid = fitQuadFree(r_ratio, kMdLo, kMdHi,
                                 kColMid, 2 /* dashed */,
                                 (tag + "_mid").c_str());

    // ---------- Scaled SIM recipe (k only) via INTEGRAL MATCHING ----------
    // k_step5_mid = I_expRatio[0.15, 0.70] / I_simShape_step5[0.15, 0.70]
    // Match k via integral over the MID region (where sim MID was originally
    // fitted). The SAME k is then applied to HIGH part of simShape_step5 —
    // HIGH becomes a sim-predicted extrapolation, NOT separately re-normalized.
    // I_expRatio uses (ratio(b) * binwidth(b)) summed over bins fully in
    //   [0.15, 0.70]; I_simShape uses the closed-form piecewise antiderivative.
    const double I_expRatio_mid =
        histIntegralXBinWidth(r_ratio, kIntLo, kIntHi);
    const double I_simShape_mid = simShapeIntegral(kIntLo, kIntHi);
    const double sim_k = (I_simShape_mid != 0.0)
        ? I_expRatio_mid / I_simShape_mid
        : 0.0;
    const double sim_k_err = 0.0;   // integral-matching: no fit error here

    printf("\nINTEGRAL-MATCHED k_step5_mid (MID-only range [%.2f, %.2f]):\n",
           kIntLo, kIntHi);
    printf("  I_expRatio  (sum ratio(b) x binwidth(b))   = %.6f\n",
           I_expRatio_mid);
    printf("  I_simShape  (closed-form piecewise int)    = %.6f\n",
           I_simShape_mid);
    printf("  k_step5_mid = I_expRatio / I_simShape       = %.4f\n", sim_k);
    printf("  -> k_step5_mid is applied to BOTH MID and HIGH parts of "
           "simShape_step5;\n");
    printf("     HIGH becomes a sim-predicted extrapolation under this "
           "normalization.\n");

    // Build the two display TF1s using the integral-matched k value.
    TF1* f_simref_MID = new TF1(("f_simref_MID_" + tag).c_str(),
                                simMidDisplay, 0.15, SIM_JOIN_X, 1);
    f_simref_MID->SetParameter(0, sim_k);
    f_simref_MID->SetLineColor(kColSimMid);
    f_simref_MID->SetLineStyle(1);   // SOLID
    f_simref_MID->SetLineWidth(3);

    TF1* f_simref_HIGH = new TF1(("f_simref_HIGH_" + tag).c_str(),
                                 simHighDisplay, SIM_JOIN_X, kXmax, 1);
    f_simref_HIGH->SetParameter(0, sim_k);
    f_simref_HIGH->SetLineColor(kColSimHigh);
    f_simref_HIGH->SetLineStyle(1);  // SOLID
    f_simref_HIGH->SetLineWidth(3);

    // ---------- Diagnostic yields ----------
    const double y_num_full = h_num->Integral();
    const double y_den_full = h_den->Integral();
    const double y_num_hi   = sumInRange(h_num, 1.0, 1.4);
    const double y_den_hi   = sumInRange(h_den, 1.0, 1.4);
    const double y_pt2_disp_hi  = sumInRange(h_disp_PT2, 1.0, 1.4) * kTrigCorr;
    const double y_pt3_disp_hi  = sumInRange(h_disp_PT3, 1.0, 1.4);

    printf("Summed yields in [1.0, 1.4]:\n");
    printf("  PT3 (raw, %s)            : %.3f\n",
           variant_label.c_str(), y_den_hi);
    printf("  63 * PT2 (%s)            : %.3f\n",
           variant_label.c_str(), 63.0 * y_num_hi);
    printf("Display-panel summed yields in [1.0, 1.4]:\n");
    printf("  PT3 display total            : %.3f\n", y_pt3_disp_hi);
    printf("  63 * PT2 display total       : %.3f\n", y_pt2_disp_hi);
    printf("Full-range yields:\n");
    printf("  PT3 (entering ratio)         : %.3f\n", y_den_full);
    printf("  PT2 (entering ratio, raw)    : %.3f\n", y_num_full);

    printf("\nNative MID FREE quad [%.2f, %.2f]:\n", kMdLo, kMdHi);
    printf("  a = %+.4f +- %.4f\n", mid.a, mid.a_err);
    printf("  b = %+.4f +- %.4f\n", mid.b, mid.b_err);
    printf("  c = %+.4f +- %.4f\n", mid.c, mid.c_err);
    printf("  chi2/ndf = %.2f / %d = %.3f\n",
           mid.chi2_ndf * mid.ndf, mid.ndf, mid.chi2_ndf);

    // ---------- DIAGNOSTIC: exp_ratio vs scaled_sim (piecewise) ------------
    // For every bin whose CENTER falls in [0.15, 0.80], dump a table with
    // (bin, x_center, exp_ratio, exp_err, scaled_sim, exp-sim, ABOVE/BELOW).
    // scaled_sim uses sim_k (JUST-FITTED) and simShape(x_center), which is
    // piecewise: MID polynomial for x<SIM_JOIN_X (0.70), HIGH for x>=0.70.
    {
        const double diag_lo = 0.15, diag_hi = 0.80;
        printf("\n==================== DIAG: exp_ratio vs scaled_sim "
               "====================\n");
        printf("  bin   x_center   exp_ratio   exp_err   scaled_sim   "
               "exp - sim   above_or_below\n");
        printf("  ---   --------   ---------   -------   ----------   "
               "---------   --------------\n");

        int    n_above = 0,  n_below = 0;
        double sum_above = 0.0, sum_below = 0.0;

        // Track first 5 bins in HIGH portion [0.70, 0.80] for separate dump.
        struct HighRow {
            int    bin;
            double xc, er, ee, ss, diff;
            const char* tag;
        };
        std::vector<HighRow> high_rows;

        for (int b = 1; b <= r_ratio->GetNbinsX(); ++b) {
            const double xc = r_ratio->GetXaxis()->GetBinCenter(b);
            if (xc < diag_lo - 1e-9 || xc > diag_hi + 1e-9) continue;
            const double bw = r_ratio->GetBinWidth(b);
            const double er = r_ratio->GetBinContent(b);
            const double ee = r_ratio->GetBinError  (b);
            const double ss = sim_k * simShape(xc);
            const double diff = er - ss;
            const char* tag = (diff > 0.0) ? "ABOVE" : "BELOW";
            printf("  %3d   %8.4f   %9.4f   %7.4f   %10.4f   %+9.4f   %s\n",
                   b, xc, er, ee, ss, diff, tag);
            if (diff > 0.0) {
                ++n_above;
                sum_above += diff * bw;
            } else if (diff < 0.0) {
                ++n_below;
                sum_below += diff * bw;
            }
            if (xc >= SIM_JOIN_X - 1e-9 && high_rows.size() < 5) {
                high_rows.push_back({b, xc, er, ee, ss, diff, tag});
            }
        }

        const double net = sum_above + sum_below;
        printf("==================== DIAG SUMMARY ====================\n");
        printf("  N_above   = %d\n", n_above);
        printf("  N_below   = %d\n", n_below);
        printf("  sum_above = %+.6f  (positive, sum of (exp-sim)*binwidth for "
               "ABOVE bins)\n", sum_above);
        printf("  sum_below = %+.6f  (negative, sum of (exp-sim)*binwidth for "
               "BELOW bins)\n", sum_below);
        printf("  net       = sum_above + sum_below = %+.6f\n", net);
        printf("  (Normalization range is [%.2f, %.2f]; diag range is "
               "[%.2f, %.2f];\n",
               kIntLo, kIntHi, diag_lo, diag_hi);
        printf("   so net should equal the contribution from [%.2f, %.2f].)\n",
               kIntHi, diag_hi);

        printf("\n  ---- FIRST 5 bins in HIGH region [0.70, 0.80] ----\n");
        printf("  bin   x_center   exp_ratio   exp_err   scaled_sim   "
               "exp - sim   above_or_below\n");
        for (const auto& r : high_rows) {
            printf("  %3d   %8.4f   %9.4f   %7.4f   %10.4f   %+9.4f   %s\n",
                   r.bin, r.xc, r.er, r.ee, r.ss, r.diff, r.tag);
        }
        printf("======================================================\n");
    }

    // ---------- Build LEFT panel display spectra ----------
    TH1D* h_PT3_w = static_cast<TH1D*>(
        h_disp_PT3->Clone(("h_PT3_w_" + tag).c_str()));
    h_PT3_w->SetDirectory(nullptr);
    TH1D* h_PT2_w = static_cast<TH1D*>(
        h_disp_PT2->Clone(("h_PT2_w_x63_" + tag).c_str()));
    h_PT2_w->SetDirectory(nullptr);
    h_PT2_w->Scale(kTrigCorr);
    normalizeByBinWidth(h_PT3_w);
    normalizeByBinWidth(h_PT2_w);
    styleHist(h_PT3_w, kColPT3, 20);
    styleHist(h_PT2_w, kColPT2, 21);

    // ---------- Hypothetical 63 * PT2 = scaledSimFit_step5(x) * PT3(x) ----------
    TH1D* h_PT2_hypo_w = static_cast<TH1D*>(
        h_disp_PT3->Clone(("h_PT2_hypo_w_" + tag).c_str()));
    h_PT2_hypo_w->SetDirectory(nullptr);
    h_PT2_hypo_w->Reset();
    for (int b = 1; b <= h_disp_PT3->GetNbinsX(); ++b) {
        const double lo = h_disp_PT3->GetXaxis()->GetBinLowEdge(b);
        if (lo < 0.15 - 1e-9) continue;   // skip LOW region
        const double n3  = h_disp_PT3->GetBinContent(b);
        const double e3  = h_disp_PT3->GetBinError  (b);
        if (n3 <= 0.0 || !std::isfinite(n3)) continue;
        const double xc  = h_disp_PT3->GetXaxis()->GetBinCenter(b);
        const double sf  = sim_k * simShape(xc);   // scaledSimFit_step5(x)
        const double val = sf * n3;
        const double err = std::abs(sf) * e3;
        h_PT2_hypo_w->SetBinContent(b, val);
        h_PT2_hypo_w->SetBinError  (b, err);
    }
    normalizeByBinWidth(h_PT2_hypo_w);
    h_PT2_hypo_w->SetMarkerStyle(25);          // open square
    h_PT2_hypo_w->SetMarkerSize(1.2);
    h_PT2_hypo_w->SetMarkerColor(kBlue + 1);
    h_PT2_hypo_w->SetLineColor  (kBlue + 1);
    h_PT2_hypo_w->SetLineWidth  (1);

    // ---------- Canvas: 1×2 ----------
    TCanvas* c = new TCanvas(("c_" + tag).c_str(), variant_label.c_str(),
                             1700, 600);
    c->Divide(2, 1, 0.001, 0.001);

    // ----- Left: spectra -----
    c->cd(1);
    gPad->SetLogy(true);
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);

    double y_max_pt3  = h_PT3_w->GetMaximum();
    double y_max_pt2  = h_PT2_w->GetMaximum();
    double y_max_glob = std::max(y_max_pt3, y_max_pt2);
    double y_min_pos  = 1e300;
    for (int b = 1; b <= h_PT3_w->GetNbinsX(); ++b) {
        const double v3 = h_PT3_w->GetBinContent(b);
        if (v3 > 0.0 && v3 < y_min_pos) y_min_pos = v3;
        const double v2 = h_PT2_w->GetBinContent(b);
        if (v2 > 0.0 && v2 < y_min_pos) y_min_pos = v2;
    }
    if (!std::isfinite(y_min_pos) || y_min_pos > 1e299) y_min_pos = 0.1;

    const double y_max = 1.5  * y_max_glob;
    const double y_min = std::max(0.1, 0.3 * y_min_pos);

    h_PT3_w->SetTitle(
        Form("%s;M_{e^{+}e^{-}} [GeV/c^{2}];weighted entries / GeV/c^{2}",
             left_title));
    h_PT3_w->GetXaxis()->SetRangeUser(kXmin, kXmax);
    h_PT3_w->GetYaxis()->SetRangeUser(y_min, y_max);
    h_PT3_w->Draw("E1");
    h_PT2_w->Draw("E1 SAME");
    h_PT2_hypo_w->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.40, 0.68, 0.95, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.028);
    leg->AddEntry(h_PT3_w,
        Form("PT3 (#sum = %.0f)", h_disp_PT3->Integral()), "lpe");
    leg->AddEntry(h_PT2_w,
        Form("63 #upoint PT2 (red dots) - measured (#sum_{PT2} = %.0f)",
             h_disp_PT2->Integral()), "lpe");
    leg->AddEntry(h_PT2_hypo_w,
        "63 #upoint PT2 hypothetical (blue open sq., x #geq 0.15) "
        "= PT3 #times scaled SIM ratio (step5)", "lpe");
    leg->Draw();

    // ----- Right: ratio + fits -----
    c->cd(2);
    gPad->SetMargin(0.13, 0.04, 0.13, 0.10);
    r_ratio->SetTitle(
        Form("%s;M_{e^{+}e^{-}} [GeV/c^{2}];%s", right_title, right_ylabel));
    r_ratio->GetXaxis()->SetRangeUser(kXmin, kXmax);
    r_ratio->GetYaxis()->SetRangeUser(0.5, 4.0);

    r_ratio->Draw("E1");

    TLine* lref = new TLine(kXmin, 1.0, kXmax, 1.0);
    lref->SetLineStyle(2); lref->SetLineColor(kGray + 2);
    lref->Draw();

    r_ratio->Draw("E1 SAME");

    // SCALED sim ratio overlay in LOW [0.0, 0.15] — small empty circles.
    TH1D* h_sim_low_scaled = buildSimLowRatio(tag);
    if (h_sim_low_scaled) {
        h_sim_low_scaled->Scale(sim_k);   // scale by k_step5
        h_sim_low_scaled->SetMarkerStyle(24);
        h_sim_low_scaled->SetMarkerSize(1.0);
        h_sim_low_scaled->SetMarkerColor(kRed);
        h_sim_low_scaled->SetLineColor  (kRed);
        h_sim_low_scaled->SetLineWidth  (1);

        printf("[sim overlay] scaled sim ratio sample values (k_step5 = %.3f):\n",
               sim_k);
        const double sample_x[] = {0.015, 0.050, 0.100, 0.135};
        for (double xv : sample_x) {
            const int sb = h_sim_low_scaled->FindBin(xv);
            const double sv = h_sim_low_scaled->GetBinContent(sb);
            const double se = h_sim_low_scaled->GetBinError  (sb);
            const double lo = h_sim_low_scaled->GetXaxis()->GetBinLowEdge(sb);
            const double hi = h_sim_low_scaled->GetXaxis()->GetBinUpEdge (sb);
            printf("    x=%.3f  bin[%.4f, %.4f]  k*ratio = %.3f +- %.3f\n",
                   xv, lo, hi, sv, se);
        }
    }

    // Fits drawn ON TOP of markers, in the requested order
    if (mid.f)        mid.f       ->Draw("SAME");
    if (f_simref_MID) f_simref_MID->Draw("SAME");
    if (f_simref_HIGH)f_simref_HIGH->Draw("SAME");

    if (h_sim_low_scaled) h_sim_low_scaled->Draw("P SAME");

    // ---------- LOW magenta: INTEGRAL-normalized in [0, 0.05] ----------
    // Same construction as the step4 version: build a scaling factor
    //   s = I_expLOW / I_simLOW
    // over [0, 0.05] using the already-k_step5-scaled red sim circles, then
    //   k_magenta = s * k_step5.
    const double k_orig = sim_k;
    double I_simLOW = 0.0;
    if (h_sim_low_scaled) {
        for (int b = 1; b <= h_sim_low_scaled->GetNbinsX(); ++b) {
            const double xc = h_sim_low_scaled->GetXaxis()->GetBinCenter(b);
            if (xc < 0.0 || xc > 0.05 + 1e-9) continue;
            const double w = h_sim_low_scaled->GetBinWidth(b);
            const double v = h_sim_low_scaled->GetBinContent(b);
            I_simLOW += v * w;
        }
    }
    double I_expLOW = 0.0;
    for (int b = 1; b <= r_ratio->GetNbinsX(); ++b) {
        const double lo = r_ratio->GetXaxis()->GetBinLowEdge(b);
        const double hi = r_ratio->GetXaxis()->GetBinUpEdge (b);
        if (lo < 0.0 - 1e-9) continue;
        if (hi > 0.05 + 1e-9) continue;
        const double w = r_ratio->GetBinWidth(b);
        const double v = r_ratio->GetBinContent(b);
        I_expLOW += v * w;
    }
    const double factor_s   = (I_simLOW != 0.0) ? I_expLOW / I_simLOW : 0.0;
    const double k_magenta  = factor_s * k_orig;

    printf("\nLOW INTEGRAL-NORMALIZED MAGENTA (range [0.00, 0.05]):\n");
    printf("  k_step5_mid (norm: integral [%.2f, %.2f] / MID-only) = %.3f\n",
           kIntLo, kIntHi, k_orig);
    printf("  I_simLOW (red circles x binwidth, [0, 0.05])     = %.4f\n",
           I_simLOW);
    printf("  I_expLOW (exp signal ratio x binwidth, [0, 0.05])= %.4f\n",
           I_expLOW);
    printf("  factor s = I_expLOW / I_simLOW                   = %.4f\n",
           factor_s);
    printf("  k_magenta = s * k_step5                          = %.4f\n",
           k_magenta);
    if (factor_s > 0.0 && factor_s < 1.0) {
        printf("  -> magenta sits %.3f x LOWER than red sim circles in [0, 0.05]\n",
               1.0 / factor_s);
    } else if (factor_s > 1.0) {
        printf("  -> magenta sits %.3f x HIGHER than red sim circles in [0, 0.05]\n",
               factor_s);
    } else {
        printf("  -> magenta equals red sim circles in [0, 0.05] (s = 1)\n");
    }

    // Magenta line.
    TF1* f_magenta = new TF1(("f_magenta_" + tag).c_str(),
                             simShapeScaled, 0.15, 1.40, 1);
    f_magenta->SetParameter(0, k_magenta);
    f_magenta->SetLineColor(kMagenta + 1);
    f_magenta->SetLineStyle(2);
    f_magenta->SetLineWidth(3);

    f_magenta->Draw("SAME");

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.026);

    tx.SetTextColor(kColMid);
    tx.DrawLatex(0.16, 0.88,
        Form("[%.2f-%.2f] native quad FREE: a=%.3f b=%.3f c=%.3f  "
             "#chi^{2}/ndf = %.1f / %d",
             kMdLo, kMdHi,
             mid.a, mid.b, mid.c,
             mid.chi2_ndf * mid.ndf, mid.ndf));

    tx.SetTextColor(kColSimMid);
    tx.DrawLatex(0.16, 0.84,
        Form("k_{step5} (norm: integral [%.2f, %.2f] / MID-only) = %.3f "
             "(shape from sim step5 same_vertex REC; HIGH = extrapolation)",
             kIntLo, kIntHi, sim_k));

    tx.SetTextColor(kColSimHigh);
    tx.DrawLatex(0.16, 0.80,
        Form("I_{exp}^{ratio} = %.4f,  I_{sim}^{shape} = %.4f  "
             "(no chi^{2})", I_expRatio_mid, I_simShape_mid));

    if (h_sim_low_scaled) {
        tx.SetTextColor(kRed);
        tx.DrawLatex(0.16, 0.76,
            Form("scaled sim points [0.0-0.15]: "
                 "k = %.2f #upoint (sim PT2/PT3, step5 same_vertex REC)",
                 sim_k));
    }

    tx.SetTextColor(kMagenta + 1);
    tx.DrawLatex(0.16, 0.72,
        Form("k_{step5} = %.2f, s = I_{exp}^{LOW}/I_{sim}^{LOW} = %.3f, "
             "k_{magenta} = s #upoint k_{step5} = %.3f",
             k_orig, factor_s, k_magenta));
    tx.DrawLatex(0.16, 0.68,
        Form("I_{simLOW} = %.3f,  I_{expLOW} = %.3f  [0, 0.05]",
             I_simLOW, I_expLOW));
    if (factor_s > 0.0 && factor_s < 1.0) {
        tx.DrawLatex(0.16, 0.64,
            Form("#rightarrow magenta sits #times %.2f LOWER than red sim "
                 "in [0, 0.05]", 1.0 / factor_s));
    } else if (factor_s > 1.0) {
        tx.DrawLatex(0.16, 0.64,
            Form("#rightarrow magenta sits #times %.2f HIGHER than red sim "
                 "in [0, 0.05]", factor_s));
    }

    TLegend* leg_r = new TLegend(0.55, 0.50, 0.95, 0.78);
    leg_r->SetBorderSize(0); leg_r->SetFillStyle(0); leg_r->SetTextSize(0.026);
    leg_r->AddEntry(mid.f,         "native MID free quad [0.10-0.80]", "l");
    leg_r->AddEntry(f_simref_MID,
        "scaled SIM step5 (MID portion)",         "l");
    leg_r->AddEntry(f_simref_HIGH,
        "scaled SIM step5 (HIGH portion)",        "l");
    leg_r->AddEntry(f_magenta,
        "INTEGRAL-matched LOW (k_{magenta}) [0.0, 1.4]", "l");
    if (h_sim_low_scaled) {
        leg_r->AddEntry(h_sim_low_scaled,
            "scaled sim PT2/PT3 (step5 same_vertex REC, LOW)", "pe");
    }
    leg_r->Draw();

    // ---------- High-mass bin comparison ----------
    printf("\nHigh-mass comparison (raw bin contents, NOT width-normalized):\n");
    printf("    x_target   bin[lo, hi]            PT3        63*PT2_actual   63*PT2_hypo  hypo/actual\n");
    const double sample_xs[] = {1.03, 1.06, 1.10, 1.20};
    for (double xs : sample_xs) {
        const int b = h_disp_PT3->FindBin(xs);
        const double lo = h_disp_PT3->GetXaxis()->GetBinLowEdge(b);
        const double hi = h_disp_PT3->GetXaxis()->GetBinUpEdge (b);
        const double xc = h_disp_PT3->GetXaxis()->GetBinCenter (b);
        const double n3 = h_disp_PT3->GetBinContent(b);
        const double n2 = h_disp_PT2->GetBinContent(b) * kTrigCorr;
        const double sf = sim_k * simShape(xc);
        const double n2_hypo = sf * n3;
        const double rr =
            (std::abs(n2) > 1e-12) ? n2_hypo / n2 : 0.0;
        printf("    x=%.2f   [%.4f, %.4f]   %8.3f    %10.3f     %10.3f    %.3f\n",
               xs, lo, hi, n3, n2, n2_hypo, rr);
    }

    // ---------- Save ----------
    c->SaveAs((base_path + ".pdf").c_str());
    c->SaveAs((base_path + ".png").c_str());
    std::cout << "Saved: " << base_path << ".{pdf,png}\n";
}

}  // anonymous namespace

void mass_ee_pt_ratio_signal_varbins_exp_step5() {
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
    std::cout << "  epem : " << kFileEm << " (" << t_em->GetEntries() << " entries)\n";
    std::cout << "  epep : " << kFilePp << " (" << t_pp->GetEntries() << " entries)\n";
    std::cout << "  emem : " << kFileMm << " (" << t_mm->GetEntries() << " entries)\n";
    std::cout << "Binning (variable, " << kNb << " bins):\n";
    std::cout << "  [0.00, 0.05] : 20 x 2.5 MeV (DENSE)\n";
    std::cout << "  [0.05, 1.40] : 38 x ~35.526 MeV (UNIFORM)\n";
    std::cout << "Trigger cuts: PT3='" << kCutPT3 << "'  PT2='" << kCutPT2 << "'\n";
    std::cout << "Correction factor (PT2 downscale 64 -> factor): "
              << kTrigCorr << "\n";
    std::cout << "SIM (step5 same_vertex) constants:\n";
    std::cout << "  SIM_JOIN_X = " << SIM_JOIN_X
              << "  SIM_HI = (" << SIM_HI_A << ", " << SIM_HI_B
              << ", " << SIM_HI_C << ")"
              << "  SIM_MID = (" << SIM_MID_B << ", " << SIM_MID_C << ")\n";
    std::cout << "  SIM_JOIN_Y = HIGH(0.70) = " << SIM_JOIN_Y << "\n\n";

    gSystem->mkdir("plots/output", true);

    // ---------- Build raw mass histograms per trigger ----------
    TH1D* h_em_PT3 = drawWithCut(t_em, "m_ee", kCutPT3, "h_em_PT3");
    TH1D* h_em_PT2 = drawWithCut(t_em, "m_ee", kCutPT2, "h_em_PT2");
    TH1D* h_pp_PT3 = drawWithCut(t_pp, "m_ee", kCutPT3, "h_pp_PT3");
    TH1D* h_pp_PT2 = drawWithCut(t_pp, "m_ee", kCutPT2, "h_pp_PT2");
    TH1D* h_mm_PT3 = drawWithCut(t_mm, "m_ee", kCutPT3, "h_mm_PT3");
    TH1D* h_mm_PT2 = drawWithCut(t_mm, "m_ee", kCutPT2, "h_mm_PT2");

    std::cout << "Raw yields:\n";
    std::cout << "  PT3: epem=" << h_em_PT3->Integral()
              << "  ++="       << h_pp_PT3->Integral()
              << "  --="       << h_mm_PT3->Integral() << "\n";
    std::cout << "  PT2: epem=" << h_em_PT2->Integral()
              << "  ++="       << h_pp_PT2->Integral()
              << "  --="       << h_mm_PT2->Integral() << "\n";

    // ---------- CB per trigger, signal per trigger ----------
    TH1D* h_cb_PT3  = makeCB(h_pp_PT3, h_mm_PT3, "h_cb_PT3");
    TH1D* h_cb_PT2  = makeCB(h_pp_PT2, h_mm_PT2, "h_cb_PT2");
    TH1D* h_sig_PT3 = makeSignal(h_em_PT3, h_cb_PT3, "h_sig_PT3");
    TH1D* h_sig_PT2 = makeSignal(h_em_PT2, h_cb_PT2, "h_sig_PT2");

    std::cout << "Signal yields (epem - CB):\n";
    std::cout << "  PT3 signal = " << h_sig_PT3->Integral() << "\n";
    std::cout << "  PT2 signal = " << h_sig_PT2->Integral() << "\n";

    // =====================================================================
    // OUTPUT: SIGNAL (CB-subtracted) — single output, single PDF/PNG
    // =====================================================================
    makeOnePlot(
        h_sig_PT2, h_sig_PT3,
        h_sig_PT3, h_sig_PT2,
        "m_{ee} CB-subtracted signal (PT3 black; 63 #upoint PT2 red)",
        "63 #upoint N_{PT2}^{sig} / N_{PT3}^{sig} (exp, step5)",
        "63 #upoint N_{PT2}^{sig} / N_{PT3}^{sig}",
        "plots/output/mass_ee_pt_ratio_varbins_exp_signal_step5",
        "signal_step5",
        "SIGNAL (CB-subtracted, step5 same_vertex)");

    f_em->Close();
    f_pp->Close();
    f_mm->Close();
}
