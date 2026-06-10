// _explore_high_binning.C
// =========================================================================
// Exploratory macro: scans candidate variable-width binning configurations
// for the HIGH region [0.80, 1.40] and the LOW-fade [0.10, 0.14] region,
// computes per-bin 63 * sig_PT2 / sig_PT3 ratios with full CB-aware error
// propagation, and ranks them by smoothness around the MID pol0 ~1.82.
//
// Usage: root -l -b -q plots/_explore_high_binning.C
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
constexpr double kMidRef   = 1.82;   // MID pol0 ratio target

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

// Run one configuration: print per-bin table; return (deviation, spread, bad).
struct ConfigStats {
    double deviation;
    double spread;
    int    nbad;
    int    ndef;
};

ConfigStats runConfig(TTree* t_em, TTree* t_pp, TTree* t_mm,
                      const std::string& tag,
                      const std::string& label,
                      const std::vector<double>& edges) {
    printf("\n==================== %s ====================\n", label.c_str());
    printf("Edges:");
    for (double e : edges) printf(" %.3f", e);
    printf("   (nbins=%zu)\n", edges.size() - 1);

    // Build histograms over this binning.
    TH1D* h_em_PT3 = drawWithCut(t_em, "m_ee", kCutPT3,
                                 "h_em_PT3_" + tag, edges);
    TH1D* h_em_PT2 = drawWithCut(t_em, "m_ee", kCutPT2,
                                 "h_em_PT2_" + tag, edges);
    TH1D* h_pp_PT3 = drawWithCut(t_pp, "m_ee", kCutPT3,
                                 "h_pp_PT3_" + tag, edges);
    TH1D* h_pp_PT2 = drawWithCut(t_pp, "m_ee", kCutPT2,
                                 "h_pp_PT2_" + tag, edges);
    TH1D* h_mm_PT3 = drawWithCut(t_mm, "m_ee", kCutPT3,
                                 "h_mm_PT3_" + tag, edges);
    TH1D* h_mm_PT2 = drawWithCut(t_mm, "m_ee", kCutPT2,
                                 "h_mm_PT2_" + tag, edges);

    TH1D* h_cb_PT3  = makeCB     (h_pp_PT3, h_mm_PT3, "h_cb_PT3_"  + tag);
    TH1D* h_cb_PT2  = makeCB     (h_pp_PT2, h_mm_PT2, "h_cb_PT2_"  + tag);
    TH1D* h_sig_PT3 = makeSignal (h_em_PT3, h_cb_PT3, "h_sig_PT3_" + tag);
    TH1D* h_sig_PT2 = makeSignal (h_em_PT2, h_cb_PT2, "h_sig_PT2_" + tag);

    printf("  %-6s  %-14s %12s %12s %10s %10s %s\n",
           "bin#", "[lo, hi]", "sig_PT3", "sig_PT2", "ratio", "sigma", "flag");

    double sum_dev   = 0.0;
    double sum_sig   = 0.0;
    int    ndef      = 0;
    int    nbad      = 0;

    for (int b = 1; b <= h_sig_PT3->GetNbinsX(); ++b) {
        const double lo = h_sig_PT3->GetXaxis()->GetBinLowEdge(b);
        const double hi = h_sig_PT3->GetXaxis()->GetBinUpEdge (b);
        const double n3 = h_sig_PT3->GetBinContent(b);
        const double e3 = h_sig_PT3->GetBinError  (b);
        const double n2 = h_sig_PT2->GetBinContent(b);
        const double e2 = h_sig_PT2->GetBinError  (b);

        std::string flag = "OK";
        double r = 0.0, s = 0.0;
        bool defined = false;

        if (!std::isfinite(n3) || !std::isfinite(n2) || n3 == 0.0) {
            flag = "BAD(zero/NaN den)";
            ++nbad;
        } else {
            r = kTrigCorr * n2 / n3;
            const double rel2 =
                (n2 != 0.0 ? std::pow(e2 / n2, 2) : 0.0) +
                std::pow(e3 / n3, 2);
            s = std::abs(r) * std::sqrt(rel2);
            if (r <= 0.0) {
                flag = "BAD(neg/zero ratio)";
                ++nbad;
            } else {
                defined = true;
            }
        }

        if (defined) {
            sum_dev += std::abs(r - kMidRef);
            sum_sig += s;
            ++ndef;
        }

        printf("  %-6d  [%.3f, %.3f] %12.3f %12.3f %10.3f %10.3f  %s\n",
               b, lo, hi, n3, n2, r, s, flag.c_str());
    }

    ConfigStats st;
    st.deviation = (ndef > 0) ? sum_dev / ndef : 1e9;
    st.spread    = (ndef > 0) ? sum_sig / ndef : 1e9;
    st.nbad      = nbad;
    st.ndef      = ndef;

    printf("  >>> ndef=%d  nbad=%d  mean|ratio-1.82|=%.3f  mean(sigma)=%.3f\n",
           st.ndef, st.nbad, st.deviation, st.spread);

    delete h_em_PT3; delete h_em_PT2;
    delete h_pp_PT3; delete h_pp_PT2;
    delete h_mm_PT3; delete h_mm_PT2;
    delete h_cb_PT3; delete h_cb_PT2;
    delete h_sig_PT3; delete h_sig_PT2;

    return st;
}

}  // anonymous namespace

void _explore_high_binning() {
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
    std::cout << "MID pol0 reference ratio = " << kMidRef << "\n";

    // -----------------------------------------------------------------
    // HIGH region [0.80, 1.40] candidate configurations.
    // -----------------------------------------------------------------
    struct Cfg {
        std::string             tag;
        std::string             label;
        std::vector<double>     edges;
    };
    std::vector<Cfg> high_cfgs = {
        {"A", "HIGH Config A: 2 bins [0.80, 1.10, 1.40]",
            {0.80, 1.10, 1.40}},
        {"B", "HIGH Config B: 3 bins [0.80, 0.95, 1.10, 1.40]",
            {0.80, 0.95, 1.10, 1.40}},
        {"C", "HIGH Config C: 4 bins [0.80, 0.95, 1.05, 1.20, 1.40]",
            {0.80, 0.95, 1.05, 1.20, 1.40}},
        {"D", "HIGH Config D: 5 bins [0.80, 0.90, 1.00, 1.10, 1.20, 1.40]",
            {0.80, 0.90, 1.00, 1.10, 1.20, 1.40}},
        {"E", "HIGH Config E: 3 asymmetric [0.80, 0.90, 1.10, 1.40]",
            {0.80, 0.90, 1.10, 1.40}},
        {"F", "HIGH Config F: 2 bins [0.80, 0.95, 1.40]",
            {0.80, 0.95, 1.40}},
    };

    std::vector<ConfigStats> high_stats;
    for (const auto& c : high_cfgs) {
        high_stats.push_back(runConfig(t_em, t_pp, t_mm, c.tag, c.label,
                                       c.edges));
    }

    // -----------------------------------------------------------------
    // LOW-fade [0.10, 0.14] candidates.
    // -----------------------------------------------------------------
    std::vector<Cfg> low_cfgs = {
        {"L1", "LOW-fade Option L1: 2 x 20 MeV [0.10, 0.12, 0.14]",
            {0.10, 0.12, 0.14}},
        {"L2", "LOW-fade Option L2: 1 x 40 MeV [0.10, 0.14]",
            {0.10, 0.14}},
    };

    std::vector<ConfigStats> low_stats;
    for (const auto& c : low_cfgs) {
        low_stats.push_back(runConfig(t_em, t_pp, t_mm, c.tag, c.label,
                                      c.edges));
    }

    // -----------------------------------------------------------------
    // Recommendation: HIGH — pick config with zero BAD bins, then min
    // deviation; tie-break by fewest bins.
    // -----------------------------------------------------------------
    printf("\n==================== SUMMARY ====================\n");
    printf("%-6s  %-50s %5s %5s %10s %10s\n",
           "tag", "label", "nbin", "nbad", "dev", "spread");
    for (size_t i = 0; i < high_cfgs.size(); ++i) {
        printf("%-6s  %-50s %5zu %5d %10.3f %10.3f\n",
               high_cfgs[i].tag.c_str(),
               high_cfgs[i].label.c_str(),
               high_cfgs[i].edges.size() - 1,
               high_stats[i].nbad,
               high_stats[i].deviation,
               high_stats[i].spread);
    }
    for (size_t i = 0; i < low_cfgs.size(); ++i) {
        printf("%-6s  %-50s %5zu %5d %10.3f %10.3f\n",
               low_cfgs[i].tag.c_str(),
               low_cfgs[i].label.c_str(),
               low_cfgs[i].edges.size() - 1,
               low_stats[i].nbad,
               low_stats[i].deviation,
               low_stats[i].spread);
    }

    int best_high = -1;
    double best_score = 1e18;
    for (size_t i = 0; i < high_cfgs.size(); ++i) {
        if (high_stats[i].nbad > 0) continue;
        // Score: deviation + 0.25 * spread, lighter penalty for more bins.
        const double score = high_stats[i].deviation
                           + 0.25 * high_stats[i].spread
                           + 0.05 * (high_cfgs[i].edges.size() - 1);
        if (score < best_score) {
            best_score = score;
            best_high  = static_cast<int>(i);
        }
    }
    int best_low = -1;
    double best_low_score = 1e18;
    for (size_t i = 0; i < low_cfgs.size(); ++i) {
        if (low_stats[i].nbad > 0) continue;
        const double score = low_stats[i].deviation
                           + 0.25 * low_stats[i].spread;
        if (score < best_low_score) {
            best_low_score = score;
            best_low = static_cast<int>(i);
        }
    }

    printf("\nRECOMMENDATION:\n");
    if (best_high >= 0) {
        printf("  HIGH  : %s   [score=%.3f]\n",
               high_cfgs[best_high].label.c_str(), best_score);
    } else {
        printf("  HIGH  : (no fully-defined configuration — all produced BAD bins)\n");
    }
    if (best_low >= 0) {
        printf("  LOW   : %s   [score=%.3f]\n",
               low_cfgs[best_low].label.c_str(), best_low_score);
    } else {
        printf("  LOW   : (no fully-defined option — all produced BAD bins)\n");
    }

    f_em->Close();
    f_pp->Close();
    f_mm->Close();
}
