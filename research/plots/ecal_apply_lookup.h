// research/plots/ecal_apply_lookup.h — bilinear / trilinear lookup of the
// ECAL energy-scale correction s(E_γ, θ_γ [, φ_local]) from a TH2D / TH3D
// map, plus the post-hoc mass correction formula used to compute m_corrected
// from per-event ntuple variables.
//
//   m²_corr = m_ee² + 2·s·E_γ·D                  (gamma_D form)
//   m²_corr = (1 − s)·m_ee² + s·m²               (closed form, no D needed)
//
// The class supports two map shapes:
//   2D  h_s_<fl>     — TH2D over (E_γ, θ_γ)
//   3D  h_s_<fl>_3d  — TH3D over (E_γ, θ_γ, φ_local)
// loadAuto() detects which is available. Lookup is bilinear / trilinear,
// with empty (zero-content) corners replaced by the mean of non-empty
// neighbours so sparse cells don't poison the surrounding region.
//
// φ convention (3D map): full φ ∈ [0°, 360°) — per-sector calibration.
// Each of HADES' 6 ECAL sectors is allowed its own correction map row
// (sectors typically have different gain/quality, so folding mod 60°
// would average out real per-sector differences).

#ifndef ECAL_APPLY_LOOKUP_H
#define ECAL_APPLY_LOOKUP_H

#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TString.h>
#include <cmath>
#include <iostream>
#include <string>

class EcalLookup {
public:
    EcalLookup() = default;

    // Load a 2D map h_s_<fl> from a ROOT file. Returns false on failure.
    bool load(const std::string& fpath, const std::string& hist_name) {
        return openFile(fpath) && grab2D(hist_name);
    }

    // Load a 3D map h_s_<fl>_3d. Returns false on failure.
    bool load3D(const std::string& fpath, const std::string& hist_name) {
        return openFile(fpath) && grab3D(hist_name);
    }

    // Try to open as 3D first, fall back to 2D. Used for callers that
    // accept either dimensionality.
    bool loadAuto(const std::string& fpath,
                  const std::string& name_2d,
                  const std::string& name_3d) {
        if (!openFile(fpath)) return false;
        if (grab3D(name_3d)) return true;
        return grab2D(name_2d);
    }

    // Normalize φ to [0, 360°). Full-range mapping per sector.
    static double phiNorm(double phi_deg) {
        const double mod = std::fmod(phi_deg, 360.0);
        return mod < 0.0 ? mod + 360.0 : mod;
    }

    // Look up s. Uses 3D if loaded (and phi_deg given), otherwise 2D.
    double s(double E_GeV, double theta_deg, double phi_deg = -999.0) const {
        if (h3_) {
            const double phi_n = (phi_deg > -900.0)
                ? phiNorm(phi_deg) : 180.0;   // arbitrary fallback if unset
            return s3D(E_GeV, theta_deg, phi_n);
        }
        if (h2_) return s2D(E_GeV, theta_deg);
        return 1.0;
    }

    // Closed-form correction — does NOT require gamma_D. Identity
    //   m²_corr = (1 − s) · m_ee² + s · m_uncorr²
    // is equivalent to m²_corr = m_ee² + 2·s·E_γ·D when the uncorrected
    // m² satisfies m² = m_ee² + 2·E_γ·D.
    double m_corrected_closed(double m_ee_GeV, double m_uncorr_GeV,
                              double E_GeV, double theta_deg,
                              double phi_deg = -999.0) const {
        const double scale = s(E_GeV, theta_deg, phi_deg);
        const double m_ee2 = m_ee_GeV * m_ee_GeV;
        const double m2    = m_uncorr_GeV * m_uncorr_GeV;
        const double m2c   = (1.0 - scale) * m_ee2 + scale * m2;
        return (m2c > 0.0) ? std::sqrt(m2c) : 0.0;
    }

    // Apply the post-hoc correction (gamma_D form).
    double m_corrected(double m_ee_GeV, double E_GeV, double D_GeV,
                       double theta_deg, double phi_deg = -999.0) const {
        const double scale = s(E_GeV, theta_deg, phi_deg);
        const double m2 = m_ee_GeV * m_ee_GeV
                          + 2.0 * scale * E_GeV * D_GeV;
        return (m2 > 0.0) ? std::sqrt(m2) : 0.0;
    }

    bool ready()    const { return h2_ != nullptr || h3_ != nullptr; }
    bool is3D()     const { return h3_ != nullptr; }
    const TH2D* hist2D() const { return h2_; }
    const TH3D* hist3D() const { return h3_; }

private:
    bool openFile(const std::string& fpath) {
        f_ = TFile::Open(fpath.c_str(), "READ");
        if (!f_ || f_->IsZombie()) {
            std::cerr << "EcalLookup: cannot open " << fpath << "\n";
            return false;
        }
        return true;
    }
    bool grab2D(const std::string& name) {
        h2_ = (TH2D*)f_->Get(name.c_str());
        if (!h2_) {
            std::cerr << "EcalLookup: TH2D '" << name << "' not found\n";
            return false;
        }
        h2_->SetDirectory(nullptr);
        return true;
    }
    bool grab3D(const std::string& name) {
        h3_ = (TH3D*)f_->Get(name.c_str());
        if (!h3_) return false;
        h3_->SetDirectory(nullptr);
        return true;
    }

    // ---- 2D bilinear with empty-corner fill --------------------------------
    double s2D(double E_GeV, double theta_deg) const {
        const TAxis* ax = h2_->GetXaxis();
        const TAxis* ay = h2_->GetYaxis();
        const double x = clampToAxis(ax, E_GeV);
        const double y = clampToAxis(ay, theta_deg);

        const int bx = ax->FindBin(x);
        const int by = ay->FindBin(y);
        const double cx = ax->GetBinCenter(bx);
        const double cy = ay->GetBinCenter(by);
        const int ix0 = (x < cx) ? bx - 1 : bx;
        const int iy0 = (y < cy) ? by - 1 : by;
        const int ix1 = ix0 + 1;
        const int iy1 = iy0 + 1;

        auto val = [&](int ix, int iy) {
            ix = std::max(1, std::min(ax->GetNbins(), ix));
            iy = std::max(1, std::min(ay->GetNbins(), iy));
            return h2_->GetBinContent(ix, iy);
        };

        double v[4] = {val(ix0,iy0), val(ix1,iy0), val(ix0,iy1), val(ix1,iy1)};
        fillEmpty(v, 4);
        if (v[0] <= 0 && v[1] <= 0 && v[2] <= 0 && v[3] <= 0) return 1.0;

        const double x0 = ax->GetBinCenter(std::max(1, std::min(ax->GetNbins(), ix0)));
        const double x1 = ax->GetBinCenter(std::max(1, std::min(ax->GetNbins(), ix1)));
        const double y0 = ay->GetBinCenter(std::max(1, std::min(ay->GetNbins(), iy0)));
        const double y1 = ay->GetBinCenter(std::max(1, std::min(ay->GetNbins(), iy1)));
        const double dx = (x1 > x0) ? (x - x0) / (x1 - x0) : 0.0;
        const double dy = (y1 > y0) ? (y - y0) / (y1 - y0) : 0.0;

        const double a = v[0] * (1 - dx) + v[1] * dx;
        const double b = v[2] * (1 - dx) + v[3] * dx;
        return a * (1 - dy) + b * dy;
    }

    // ---- 3D trilinear with empty-corner fill -------------------------------
    double s3D(double E_GeV, double theta_deg, double phi_local) const {
        const TAxis* ax = h3_->GetXaxis();
        const TAxis* ay = h3_->GetYaxis();
        const TAxis* az = h3_->GetZaxis();
        const double x = clampToAxis(ax, E_GeV);
        const double y = clampToAxis(ay, theta_deg);
        const double z = clampToAxis(az, phi_local);

        const int bx = ax->FindBin(x);
        const int by = ay->FindBin(y);
        const int bz = az->FindBin(z);
        const int ix0 = (x < ax->GetBinCenter(bx)) ? bx - 1 : bx;
        const int iy0 = (y < ay->GetBinCenter(by)) ? by - 1 : by;
        const int iz0 = (z < az->GetBinCenter(bz)) ? bz - 1 : bz;
        const int ix1 = ix0 + 1, iy1 = iy0 + 1, iz1 = iz0 + 1;

        auto val = [&](int ix, int iy, int iz) {
            ix = std::max(1, std::min(ax->GetNbins(), ix));
            iy = std::max(1, std::min(ay->GetNbins(), iy));
            iz = std::max(1, std::min(az->GetNbins(), iz));
            return h3_->GetBinContent(ix, iy, iz);
        };

        // 8 corners.
        double v[8] = {
            val(ix0, iy0, iz0), val(ix1, iy0, iz0),
            val(ix0, iy1, iz0), val(ix1, iy1, iz0),
            val(ix0, iy0, iz1), val(ix1, iy0, iz1),
            val(ix0, iy1, iz1), val(ix1, iy1, iz1)
        };
        fillEmpty(v, 8);
        bool any = false;
        for (int i = 0; i < 8; ++i) if (v[i] > 0) { any = true; break; }
        if (!any) return 1.0;

        const double x0 = ax->GetBinCenter(std::max(1, std::min(ax->GetNbins(), ix0)));
        const double x1 = ax->GetBinCenter(std::max(1, std::min(ax->GetNbins(), ix1)));
        const double y0 = ay->GetBinCenter(std::max(1, std::min(ay->GetNbins(), iy0)));
        const double y1 = ay->GetBinCenter(std::max(1, std::min(ay->GetNbins(), iy1)));
        const double z0 = az->GetBinCenter(std::max(1, std::min(az->GetNbins(), iz0)));
        const double z1 = az->GetBinCenter(std::max(1, std::min(az->GetNbins(), iz1)));
        const double dx = (x1 > x0) ? (x - x0) / (x1 - x0) : 0.0;
        const double dy = (y1 > y0) ? (y - y0) / (y1 - y0) : 0.0;
        const double dz = (z1 > z0) ? (z - z0) / (z1 - z0) : 0.0;

        // Trilinear: interpolate along x, then y, then z.
        const double c00 = v[0] * (1 - dx) + v[1] * dx;
        const double c10 = v[2] * (1 - dx) + v[3] * dx;
        const double c01 = v[4] * (1 - dx) + v[5] * dx;
        const double c11 = v[6] * (1 - dx) + v[7] * dx;
        const double c0  = c00 * (1 - dy) + c10 * dy;
        const double c1  = c01 * (1 - dy) + c11 * dy;
        return c0 * (1 - dz) + c1 * dz;
    }

    static double clampToAxis(const TAxis* a, double v) {
        const double lo = a->GetBinLowEdge(1);
        const double hi = a->GetBinUpEdge(a->GetNbins());
        // tiny inset to keep FindBin from returning overflow on the edge
        return std::max(lo, std::min(hi * (1 - 1e-9), v));
    }

    // Replace ≤0 entries with the mean of >0 entries (in-place).
    static void fillEmpty(double* v, int n) {
        double sum = 0.0;
        int    cnt = 0;
        for (int i = 0; i < n; ++i) if (v[i] > 0) { sum += v[i]; ++cnt; }
        if (cnt == 0) return;
        const double mean = sum / cnt;
        for (int i = 0; i < n; ++i) if (v[i] <= 0) v[i] = mean;
    }

    TFile* f_  = nullptr;
    TH2D*  h2_ = nullptr;
    TH3D*  h3_ = nullptr;
};

#endif
