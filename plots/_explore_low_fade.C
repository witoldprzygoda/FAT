// _explore_low_fade.C
// =========================================================================
// Exploratory macro: scans candidate variable-width binning configurations
// for the LOW-fade region [0.07, 0.14], computing per-bin
// 63 * sig_PT2 / sig_PT3 ratios with CB-aware error propagation.
//
// Each candidate combines DENSE 2.5 MeV bins in [0, 0.07] (28 bins) with a
// variable-width tail in [0.07, 0.14]. We then inspect ONLY the variable
// portion and check smoothness against:
//   - the mean of the last 2 dense bins (centers ~0.0625, 0.0675), and
//   - the MID pol0 reference ratio = 1.821.
//
// Usage: root -l -b -q plots/_explore_low_fade.C
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
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

constexpr const char* kCutPT3 = "trigbit==8192";
constexpr const char* kCutPT2 = "trigbit==4096";

constexpr double kTrigCorr = 63.0;
constexpr double kMidRef   = 1.821;   // MID pol0 ratio target
constexpr double kRatioLo  = 1.5;     // reasonable per-bin range
constexpr double kRatioHi  = 2.2;

// Build edge array: 2.5 MeV bins in [0, 0.07] (28 bins, edges 0..0.070)
// followed by the variable portion (skipping the duplicate 0.070).
std::vector<double> buildEdges(const std::vector<double>& var_edges) {
    std::vector<double> edges;
    edges.reserve(28 + var_edges.size());
    for (int i = 0; i <= 28; ++i) {
        edges.push_back(0.0025 * i);  // 0.000, 0.0025, ... 0.070
    }
    // var_edges starts at 0.07 — skip it to avoid duplicate.
    for (size_t i = 1; i < var_edges.size(); ++i) {
        edges.push_back(var_edges[i]);
    }
    return edges;
}

// Fill a TH1D over an arbitrary edge array.
TH1D* drawWithCut(TTree* t, const char* expr, const char* cut,
                  const std::string& name,
                  const std::vector<double>& edges) {
    TH1D* h = new TH1D(name.c_str(), "", static_cast<int>(edges.size()) - 1,
                      edges.data());
    h->Sumw2();
    t->Draw((std::string(expr) + ">>" + name).c_str(), cut, "goff");
    h->SetDirectory(nullptr);
    return h;
}

// CB = 2*sqrt(N_++ * N_--), per-bin, with relative-error propagation.
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
        const double err = val * 0.5 *
            std::sqrt(rel_pp * rel_pp + rel_mm * rel_mm);
        cb->SetBinContent(b, val);
        cb->SetBinError  (b, err);
    }
    return cb;
}

// signal = epem - CB with quadrature errors.
TH1D* makeSignal(TH1D* h_em, TH1D* h_cb, const std::string& name) {
    TH1D* sig = static_cast<TH1D*>(h_em->Clone(name.c_str()));
    sig->SetDirectory(nullptr);
    sig->Reset();
    for (int b = 1; b <= h_em->GetNbinsX(); ++b) {
        const double n_em = h_em->GetBinContent(b);
        const double e_em = h_em->GetBinError(b);
        const double n_cb = h_cb->GetBinContent(b);
        const double e_cb = h_cb->GetBinError(b);
        sig->SetBinContent(b, n_em - n_cb);
        sig->SetBinError  (b, std::sqrt(e_em * e_em + e_cb * e_cb));
    }
    return sig;
}

// Per-bin info we keep for each variable bin.
struct BinInfo {
    double xc;
    double width;
    double sig3;
    double sig2;
    double ratio;
    double sigma;
    bool   defined;
    bool   out_of_range;       // ratio outside [1.5, 2.2]
    bool   jump_vs_prev;       // > 2 sigma from previous variable bin
};

struct ConfigStats {
    int    nbad;
    int    n_out_of_range;
    int    n_jumps;
    double dev_dense;          // |first var ratio - dense_ref|
    double dev_mid;            // |last  var ratio - kMidRef|
    double mean_sigma;
    int    ndef;
    std::vector<BinInfo> bins;
};

ConfigStats runConfig(TTree* t_em, TTree* t_pp, TTree* t_mm,
                      const std::string& tag,
                      const std::string& label,
                      const std::vector<double>& var_edges,
                      double& dense_ref_out,
                      double& dense_ref_sigma_out) {
    printf("\n==================== %s ====================\n", label.c_str());
    printf("Variable edges in [0.07, 0.14]:");
    for (double e : var_edges) printf(" %.4f", e);
    printf("   (nvar=%zu)\n", var_edges.size() - 1);

    const std::vector<double> edges = buildEdges(var_edges);

    TH1D* h_em_PT3 = drawWithCut(t_em, "m_ee", kCutPT3,
                                 "h_em_PT3_low_" + tag, edges);
    TH1D* h_em_PT2 = drawWithCut(t_em, "m_ee", kCutPT2,
                                 "h_em_PT2_low_" + tag, edges);
    TH1D* h_pp_PT3 = drawWithCut(t_pp, "m_ee", kCutPT3,
                                 "h_pp_PT3_low_" + tag, edges);
    TH1D* h_pp_PT2 = drawWithCut(t_pp, "m_ee", kCutPT2,
                                 "h_pp_PT2_low_" + tag, edges);
    TH1D* h_mm_PT3 = drawWithCut(t_mm, "m_ee", kCutPT3,
                                 "h_mm_PT3_low_" + tag, edges);
    TH1D* h_mm_PT2 = drawWithCut(t_mm, "m_ee", kCutPT2,
                                 "h_mm_PT2_low_" + tag, edges);

    TH1D* h_cb_PT3  = makeCB     (h_pp_PT3, h_mm_PT3, "h_cb_PT3_low_"  + tag);
    TH1D* h_cb_PT2  = makeCB     (h_pp_PT2, h_mm_PT2, "h_cb_PT2_low_"  + tag);
    TH1D* h_sig_PT3 = makeSignal (h_em_PT3, h_cb_PT3, "h_sig_PT3_low_" + tag);
    TH1D* h_sig_PT2 = makeSignal (h_em_PT2, h_cb_PT2, "h_sig_PT2_low_" + tag);

    auto perBinRatio = [&](int b, double& r, double& s) -> bool {
        const double n3 = h_sig_PT3->GetBinContent(b);
        const double e3 = h_sig_PT3->GetBinError  (b);
        const double n2 = h_sig_PT2->GetBinContent(b);
        const double e2 = h_sig_PT2->GetBinError  (b);
        if (!std::isfinite(n3) || !std::isfinite(n2) || n3 == 0.0) {
            r = 0.0; s = 0.0;
            return false;
        }
        r = kTrigCorr * n2 / n3;
        const double rel2 =
            (n2 != 0.0 ? std::pow(e2 / n2, 2) : 0.0) +
            std::pow(e3 / n3, 2);
        s = std::abs(r) * std::sqrt(rel2);
        return (r > 0.0);
    };

    // Reference from the last 2 dense bins (b=27 -> [0.065,0.0675],
    // b=28 -> [0.0675,0.07]).  Inverse-sigma^2 weighted mean.
    double w_sum = 0.0, w_val = 0.0, w_var = 0.0;
    for (int b = 27; b <= 28; ++b) {
        double r, s;
        if (perBinRatio(b, r, s) && s > 0.0) {
            const double w = 1.0 / (s * s);
            w_sum += w;
            w_val += w * r;
            w_var += w * w * s * s;
        }
    }
    double dense_ref = (w_sum > 0.0) ? w_val / w_sum : 0.0;
    double dense_ref_sigma = (w_sum > 0.0) ? std::sqrt(w_var) / w_sum : 0.0;
    dense_ref_out = dense_ref;
    dense_ref_sigma_out = dense_ref_sigma;
    printf("  dense-side ref (mean of last 2 dense bins) = "
           "%.3f +/- %.3f\n", dense_ref, dense_ref_sigma);
    printf("  MID-side ref (pol0)                         = %.3f\n",
           kMidRef);

    printf("  %-4s  %-16s %8s %12s %12s %8s %8s  %s\n",
           "bin#", "[lo, hi]", "width", "sig_PT3", "sig_PT2",
           "ratio", "sigma", "flag");

    ConfigStats st;
    st.nbad = 0;
    st.n_out_of_range = 0;
    st.n_jumps = 0;
    st.dev_dense = 0.0;
    st.dev_mid = 0.0;
    st.mean_sigma = 0.0;
    st.ndef = 0;

    const int n_var = static_cast<int>(var_edges.size()) - 1;
    const int b0 = 29;          // first variable bin index
    double sum_sig = 0.0;

    double prev_r = 0.0;
    double prev_s = 0.0;
    bool   prev_ok = false;

    for (int k = 0; k < n_var; ++k) {
        const int b = b0 + k;
        const double lo = h_sig_PT3->GetXaxis()->GetBinLowEdge(b);
        const double hi = h_sig_PT3->GetXaxis()->GetBinUpEdge (b);
        const double n3 = h_sig_PT3->GetBinContent(b);
        const double n2 = h_sig_PT2->GetBinContent(b);

        BinInfo bi;
        bi.xc = 0.5 * (lo + hi);
        bi.width = hi - lo;
        bi.sig3 = n3;
        bi.sig2 = n2;
        bi.ratio = 0.0;
        bi.sigma = 0.0;
        bi.defined = false;
        bi.out_of_range = false;
        bi.jump_vs_prev = false;

        std::string flag = "OK";

        double r = 0.0, s = 0.0;
        if (!perBinRatio(b, r, s)) {
            flag = "BAD(zero/neg)";
            ++st.nbad;
        } else {
            bi.defined = true;
            bi.ratio = r;
            bi.sigma = s;
            if (r < kRatioLo || r > kRatioHi) {
                bi.out_of_range = true;
                ++st.n_out_of_range;
                flag = "OUT_OF_RANGE";
            }
            if (prev_ok) {
                const double diff = std::abs(r - prev_r);
                const double tot  = std::sqrt(s * s + prev_s * prev_s);
                if (tot > 0.0 && diff > 2.0 * tot) {
                    bi.jump_vs_prev = true;
                    ++st.n_jumps;
                    if (flag == "OK") flag = "JUMP";
                    else              flag += "+JUMP";
                }
            }
            prev_r = r; prev_s = s; prev_ok = true;
            sum_sig += s;
            ++st.ndef;
        }
        st.bins.push_back(bi);

        printf("  %-4d  [%.4f,%.4f] %8.4f %12.3f %12.3f %8.3f %8.3f  %s\n",
               k + 1, lo, hi, bi.width, n3, n2, r, s, flag.c_str());
    }

    if (st.ndef > 0) {
        st.mean_sigma = sum_sig / st.ndef;
        if (st.bins.front().defined) {
            st.dev_dense = std::abs(st.bins.front().ratio - dense_ref);
        }
        if (st.bins.back().defined) {
            st.dev_mid = std::abs(st.bins.back().ratio - kMidRef);
        }
    }

    printf("  >>> nvar=%d  ndef=%d  nbad=%d  n_out_of_range=%d  n_jumps=%d\n",
           n_var, st.ndef, st.nbad, st.n_out_of_range, st.n_jumps);
    printf("  >>> |first-dense_ref|=%.3f  |last-MIDref|=%.3f  mean(sigma)=%.3f\n",
           st.dev_dense, st.dev_mid, st.mean_sigma);

    delete h_em_PT3; delete h_em_PT2;
    delete h_pp_PT3; delete h_pp_PT2;
    delete h_mm_PT3; delete h_mm_PT2;
    delete h_cb_PT3; delete h_cb_PT2;
    delete h_sig_PT3; delete h_sig_PT2;

    return st;
}

}  // anonymous namespace

void _explore_low_fade() {
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
    std::cout << "MID pol0 reference ratio = " << kMidRef << "\n";
    std::cout << "Dense bins in [0, 0.07]: 28 x 2.5 MeV\n";

    struct Cfg {
        std::string             tag;
        std::string             label;
        std::vector<double>     var_edges;
    };
    std::vector<Cfg> cfgs = {
        {"T1", "T1: 1 bin x 70 MeV          [0.07, 0.14]",
            {0.07, 0.14}},
        {"T2", "T2: 2 bins x 35 MeV         [0.07, 0.105, 0.14]",
            {0.07, 0.105, 0.14}},
        {"T3", "T3: 2 bins (30+40 MeV)      [0.07, 0.10, 0.14]",
            {0.07, 0.10, 0.14}},
        {"T4", "T4: 3 bins (20+20+30 MeV)   [0.07, 0.09, 0.11, 0.14]",
            {0.07, 0.09, 0.11, 0.14}},
        {"T5", "T5: 3 bins x 23.3 MeV       [0.07, 0.0933, 0.1167, 0.14]",
            {0.07, 0.09333333, 0.11666667, 0.14}},
        {"T6", "T6: 4 bins (15,15,20,20)    [0.07, 0.085, 0.10, 0.12, 0.14]",
            {0.07, 0.085, 0.10, 0.12, 0.14}},
        {"T7", "T7: 4 bins x 17.5 MeV       [0.07, 0.0875, 0.105, 0.1225, 0.14]",
            {0.07, 0.0875, 0.105, 0.1225, 0.14}},
    };

    std::vector<ConfigStats> stats;
    std::vector<double> dense_refs;
    std::vector<double> dense_ref_sigmas;
    stats.reserve(cfgs.size());
    for (const auto& c : cfgs) {
        double dref = 0.0, dref_s = 0.0;
        stats.push_back(
            runConfig(t_em, t_pp, t_mm, c.tag, c.label, c.var_edges,
                      dref, dref_s));
        dense_refs.push_back(dref);
        dense_ref_sigmas.push_back(dref_s);
    }

    // -----------------------------------------------------------------
    // Summary table.
    // -----------------------------------------------------------------
    printf("\n==================== SUMMARY ====================\n");
    printf("%-4s  %-50s %5s %5s %5s %5s %8s %8s %8s\n",
           "tag", "label", "nvar", "nbad", "noor", "njmp",
           "dDense", "dMid", "<sig>");
    for (size_t i = 0; i < cfgs.size(); ++i) {
        printf("%-4s  %-50s %5zu %5d %5d %5d %8.3f %8.3f %8.3f\n",
               cfgs[i].tag.c_str(),
               cfgs[i].label.c_str(),
               cfgs[i].var_edges.size() - 1,
               stats[i].nbad,
               stats[i].n_out_of_range,
               stats[i].n_jumps,
               stats[i].dev_dense,
               stats[i].dev_mid,
               stats[i].mean_sigma);
    }
    printf("  (noor = bins with ratio outside [%.2f, %.2f]; "
           "njmp = bins jumping >2sigma vs previous variable bin)\n",
           kRatioLo, kRatioHi);

    // -----------------------------------------------------------------
    // Recommendation: zero BAD, zero out-of-range, zero jumps preferred;
    // score = dDense + dMid + 0.5*<sigma> + 0.05*nvar
    // -----------------------------------------------------------------
    int best = -1;
    double best_score = 1e18;
    for (size_t i = 0; i < cfgs.size(); ++i) {
        if (stats[i].nbad > 0) continue;
        const double penalty = 0.5 * stats[i].n_out_of_range +
                               0.5 * stats[i].n_jumps;
        const double score = stats[i].dev_dense + stats[i].dev_mid +
                             0.5 * stats[i].mean_sigma +
                             0.05 * (cfgs[i].var_edges.size() - 1) +
                             penalty;
        if (score < best_score) {
            best_score = score;
            best = static_cast<int>(i);
        }
    }

    printf("\nRECOMMENDATION:\n");
    if (best >= 0) {
        printf("  %s   [score=%.3f]\n",
               cfgs[best].label.c_str(), best_score);
        printf("  dense_ref=%.3f+/-%.3f   MID_ref=%.3f\n",
               dense_refs[best], dense_ref_sigmas[best], kMidRef);
    } else {
        printf("  (no fully-defined configuration — all produced BAD bins)\n");
    }

    f_em->Close();
    f_pp->Close();
    f_mm->Close();
}
