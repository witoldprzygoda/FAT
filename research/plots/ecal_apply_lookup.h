// research/plots/ecal_apply_lookup.h — bilinear lookup of ECAL energy-scale
// correction s(E_γ, θ_γ) from a TH2D map, plus the post-hoc mass correction
// formula used to compute m_corrected from per-event ntuple variables.
//
//   m²_corr = m_ee² + 2·s·E_γ·D
//   D = E_ee − p_ee·n̂_γ        (precomputed and stored as `gamma_D` per event)
//
// Out-of-bounds lookups clamp to the nearest valid bin (no extrapolation).
// Cells with content == 0 are treated as "no measurement" — the lookup
// falls back to the nearest valid neighbour. This keeps the corrected
// distribution reasonable in regions where the calibration scan had too
// few events to produce a good fit (e.g. θ > 45°, E_γ > 1.5 GeV).

#ifndef ECAL_APPLY_LOOKUP_H
#define ECAL_APPLY_LOOKUP_H

#include <TH2D.h>
#include <TFile.h>
#include <TString.h>
#include <cmath>
#include <iostream>
#include <string>

class EcalLookup {
public:
    EcalLookup() = default;

    // Load h_s_<fl> from a 2D scan ROOT file. Returns false on failure.
    bool load(const std::string& fpath, const std::string& hist_name) {
        f_ = TFile::Open(fpath.c_str(), "READ");
        if (!f_ || f_->IsZombie()) {
            std::cerr << "EcalLookup: cannot open " << fpath << "\n";
            return false;
        }
        h_ = (TH2D*)f_->Get(hist_name.c_str());
        if (!h_) {
            std::cerr << "EcalLookup: hist '" << hist_name
                      << "' not in " << fpath << "\n";
            return false;
        }
        h_->SetDirectory(nullptr);
        return true;
    }

    // Bilinear interpolation of s on the (E, θ) grid. Empty cells (content
    // == 0) are skipped and the nearest non-empty value is used.
    // Returns 1.0 if the map is unavailable (no correction).
    double s(double E_GeV, double theta_deg) const {
        if (!h_) return 1.0;
        const TAxis* ax = h_->GetXaxis();
        const TAxis* ay = h_->GetYaxis();
        const double xmin = ax->GetBinLowEdge(1);
        const double xmax = ax->GetBinUpEdge(ax->GetNbins());
        const double ymin = ay->GetBinLowEdge(1);
        const double ymax = ay->GetBinUpEdge(ay->GetNbins());
        const double x = std::max(xmin, std::min(xmax * (1 - 1e-9), E_GeV));
        const double y = std::max(ymin, std::min(ymax * (1 - 1e-9), theta_deg));

        // Find surrounding bin centers (i, j) and (i+1, j+1).
        const int bx = ax->FindBin(x);
        const int by = ay->FindBin(y);

        // Use the four neighbour bin centers for bilinear; fall back to
        // nearest non-empty if any neighbour is empty.
        auto val = [&](int ix, int iy) {
            ix = std::max(1, std::min(ax->GetNbins(), ix));
            iy = std::max(1, std::min(ay->GetNbins(), iy));
            return h_->GetBinContent(ix, iy);
        };

        const double cx = ax->GetBinCenter(bx);
        const double cy = ay->GetBinCenter(by);
        const int    ix0 = (x < cx) ? bx - 1 : bx;
        const int    iy0 = (y < cy) ? by - 1 : by;
        const int    ix1 = ix0 + 1;
        const int    iy1 = iy0 + 1;

        double v00 = val(ix0, iy0);
        double v10 = val(ix1, iy0);
        double v01 = val(ix0, iy1);
        double v11 = val(ix1, iy1);

        // Replace empty corners by an average of non-empty neighbours.
        const double good[4] = {v00, v10, v01, v11};
        double sum_good = 0.0;
        int    n_good = 0;
        for (double v : good) if (v > 0) { sum_good += v; ++n_good; }
        if (n_good == 0) return 1.0;
        const double mean = sum_good / n_good;
        if (v00 <= 0) v00 = mean;
        if (v10 <= 0) v10 = mean;
        if (v01 <= 0) v01 = mean;
        if (v11 <= 0) v11 = mean;

        // Bilinear weights.
        const double x0 = ax->GetBinCenter(std::max(1, std::min(ax->GetNbins(), ix0)));
        const double x1 = ax->GetBinCenter(std::max(1, std::min(ax->GetNbins(), ix1)));
        const double y0 = ay->GetBinCenter(std::max(1, std::min(ay->GetNbins(), iy0)));
        const double y1 = ay->GetBinCenter(std::max(1, std::min(ay->GetNbins(), iy1)));
        const double dx = (x1 > x0) ? (x - x0) / (x1 - x0) : 0.0;
        const double dy = (y1 > y0) ? (y - y0) / (y1 - y0) : 0.0;

        const double v0 = v00 * (1 - dx) + v10 * dx;
        const double v1 = v01 * (1 - dx) + v11 * dx;
        return v0 * (1 - dy) + v1 * dy;
    }

    // Closed-form correction — does NOT require gamma_D. Identity
    //   m²_corr = (1 − s) · m_ee² + s · m_uncorr²
    // is equivalent to m²_corr = m_ee² + 2·s·E_γ·D when the uncorrected
    // m² satisfies m² = m_ee² + 2·E_γ·D. Use this on ntuples produced
    // before gamma_D was added to meson_dalitz_nt.
    double m_corrected_closed(double m_ee_GeV, double m_uncorr_GeV,
                              double E_GeV, double theta_deg) const {
        const double scale = s(E_GeV, theta_deg);
        const double m_ee2 = m_ee_GeV * m_ee_GeV;
        const double m2    = m_uncorr_GeV * m_uncorr_GeV;
        const double m2c   = (1.0 - scale) * m_ee2 + scale * m2;
        return (m2c > 0.0) ? std::sqrt(m2c) : 0.0;
    }

    // Apply the post-hoc correction (gamma_D variant).
    double m_corrected(double m_ee_GeV, double E_GeV, double D_GeV,
                       double theta_deg) const {
        const double scale = s(E_GeV, theta_deg);
        const double m2 = m_ee_GeV * m_ee_GeV
                          + 2.0 * scale * E_GeV * D_GeV;
        return (m2 > 0.0) ? std::sqrt(m2) : 0.0;
    }

    bool ready() const { return h_ != nullptr; }
    const TH2D* hist() const { return h_; }

private:
    TFile* f_ = nullptr;
    TH2D*  h_ = nullptr;
};

#endif
