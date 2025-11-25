/**
 * @file HistogramBuilder.h
 * @brief Fluent interface builder for creating histograms
 *
 * This file provides a builder pattern for histogram creation with a fluent interface.
 * It allows for expressive, readable histogram creation with optional parameters.
 *
 * Example usage:
 *   auto hist = HistogramBuilder()
 *                   .name("h_theta")
 *                   .title("Theta distribution")
 *                   .bins(100, 0, 180)
 *                   .folder("angular")
 *                   .description("Scattering angle in lab frame")
 *                   .tag("proton")
 *                   .tag("angular")
 *                   .build1D();
 *
 * Benefits:
 * - Named parameters (no need to remember parameter order)
 * - Optional parameters with sensible defaults
 * - Readable, self-documenting code
 * - Compile-time type safety
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef HISTOGRAM_BUILDER_H
#define HISTOGRAM_BUILDER_H

#include <string>
#include <vector>
#include <memory>
#include <stdexcept>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include "histogram_registry.h"

// ============================================================================
// HistogramBuilder: Fluent interface for histogram creation
// ============================================================================

class HistogramBuilder {
public:

    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    HistogramBuilder() = default;

    // ------------------------------------------------------------------------
    // Basic properties
    // ------------------------------------------------------------------------

    /**
     * @brief Set histogram name (required)
     */
    HistogramBuilder& name(const std::string& n) {
        name_ = n;
        return *this;
    }

    /**
     * @brief Set histogram title (optional, defaults to name)
     */
    HistogramBuilder& title(const std::string& t) {
        title_ = t;
        has_title_ = true;
        return *this;
    }

    // ------------------------------------------------------------------------
    // Binning (1D)
    // ------------------------------------------------------------------------

    /**
     * @brief Set binning for X-axis
     */
    HistogramBuilder& bins(int n, double low, double up) {
        nbinsx_ = n;
        xlow_ = low;
        xup_ = up;
        has_binning_ = true;
        return *this;
    }

    /**
     * @brief Set binning for X-axis (alternative name)
     */
    HistogramBuilder& binsX(int n, double low, double up) {
        return bins(n, low, up);
    }

    // ------------------------------------------------------------------------
    // Binning (2D)
    // ------------------------------------------------------------------------

    /**
     * @brief Set binning for Y-axis (for 2D histograms)
     */
    HistogramBuilder& binsY(int n, double low, double up) {
        nbinsy_ = n;
        ylow_ = low;
        yup_ = up;
        has_binning_y_ = true;
        return *this;
    }

    // ------------------------------------------------------------------------
    // Binning (3D)
    // ------------------------------------------------------------------------

    /**
     * @brief Set binning for Z-axis (for 3D histograms)
     */
    HistogramBuilder& binsZ(int n, double low, double up) {
        nbinsz_ = n;
        zlow_ = low;
        zup_ = up;
        has_binning_z_ = true;
        return *this;
    }

    // ------------------------------------------------------------------------
    // Metadata
    // ------------------------------------------------------------------------

    /**
     * @brief Set ROOT folder path
     */
    HistogramBuilder& folder(const std::string& f) {
        folder_ = f;
        return *this;
    }

    /**
     * @brief Set description
     */
    HistogramBuilder& description(const std::string& d) {
        description_ = d;
        return *this;
    }

    /**
     * @brief Add a tag (can be called multiple times)
     */
    HistogramBuilder& tag(const std::string& t) {
        tags_.push_back(t);
        return *this;
    }

    /**
     * @brief Add multiple tags at once
     */
    HistogramBuilder& tags(const std::vector<std::string>& t) {
        tags_.insert(tags_.end(), t.begin(), t.end());
        return *this;
    }

    // ------------------------------------------------------------------------
    // Build methods
    // ------------------------------------------------------------------------

    /**
     * @brief Build a 1D histogram (TH1F)
     */
    std::unique_ptr<TH1F> build1D() {
        validate1D();

        const char* title_cstr = has_title_ ? title_.c_str() : name_.c_str();

        return std::make_unique<TH1F>(
            name_.c_str(),
            title_cstr,
            nbinsx_, xlow_, xup_
        );
    }

    /**
     * @brief Build a 2D histogram (TH2F)
     */
    std::unique_ptr<TH2F> build2D() {
        validate2D();

        const char* title_cstr = has_title_ ? title_.c_str() : name_.c_str();

        return std::make_unique<TH2F>(
            name_.c_str(),
            title_cstr,
            nbinsx_, xlow_, xup_,
            nbinsy_, ylow_, yup_
        );
    }

    /**
     * @brief Build a 3D histogram (TH3F)
     */
    std::unique_ptr<TH3F> build3D() {
        validate3D();

        const char* title_cstr = has_title_ ? title_.c_str() : name_.c_str();

        return std::make_unique<TH3F>(
            name_.c_str(),
            title_cstr,
            nbinsx_, xlow_, xup_,
            nbinsy_, ylow_, yup_,
            nbinsz_, zlow_, zup_
        );
    }

    // ------------------------------------------------------------------------
    // Build and register in one step
    // ------------------------------------------------------------------------

    /**
     * @brief Build 1D histogram and add to registry
     *
     * Example:
     *   HistogramBuilder()
     *       .name("h_theta")
     *       .bins(100, 0, 180)
     *       .folder("angular")
     *       .buildAndRegister1D(registry);
     */
    template<typename Registry>
    void buildAndRegister1D(Registry& registry) {
        auto hist = build1D();
        HistogramMetadata meta = buildMetadata();
        registry.add(std::move(hist), meta);
    }

    /**
     * @brief Build 2D histogram and add to registry
     */
    template<typename Registry>
    void buildAndRegister2D(Registry& registry) {
        auto hist = build2D();
        HistogramMetadata meta = buildMetadata();
        registry.add(std::move(hist), meta);
    }

    /**
     * @brief Build 3D histogram and add to registry
     */
    template<typename Registry>
    void buildAndRegister3D(Registry& registry) {
        auto hist = build3D();
        HistogramMetadata meta = buildMetadata();
        registry.add(std::move(hist), meta);
    }

    // ------------------------------------------------------------------------
    // Get metadata
    // ------------------------------------------------------------------------

    /**
     * @brief Build metadata object (useful for manual registration)
     */
    HistogramMetadata buildMetadata() const {
        HistogramMetadata meta(name_, folder_, description_);
        meta.tags = tags_;
        return meta;
    }

    // ------------------------------------------------------------------------
    // Reset builder (for reuse)
    // ------------------------------------------------------------------------

    /**
     * @brief Reset all parameters (allows builder reuse)
     */
    void reset() {
        name_.clear();
        title_.clear();
        folder_.clear();
        description_.clear();
        tags_.clear();

        has_title_ = false;
        has_binning_ = false;
        has_binning_y_ = false;
        has_binning_z_ = false;

        nbinsx_ = 100;
        xlow_ = 0.0;
        xup_ = 100.0;

        nbinsy_ = 100;
        ylow_ = 0.0;
        yup_ = 100.0;

        nbinsz_ = 100;
        zlow_ = 0.0;
        zup_ = 100.0;
    }

private:

    // ------------------------------------------------------------------------
    // Validation
    // ------------------------------------------------------------------------

    void validate1D() const {
        if (name_.empty()) {
            throw std::runtime_error("HistogramBuilder: name() must be called before build1D()");
        }
        if (!has_binning_) {
            throw std::runtime_error("HistogramBuilder: bins() must be called before build1D()");
        }
    }

    void validate2D() const {
        validate1D();
        if (!has_binning_y_) {
            throw std::runtime_error("HistogramBuilder: binsY() must be called before build2D()");
        }
    }

    void validate3D() const {
        validate2D();
        if (!has_binning_z_) {
            throw std::runtime_error("HistogramBuilder: binsZ() must be called before build3D()");
        }
    }

    // ------------------------------------------------------------------------
    // Member variables
    // ------------------------------------------------------------------------

    // Basic properties
    std::string name_;
    std::string title_;
    bool has_title_ = false;

    // Binning (X)
    int nbinsx_ = 100;
    double xlow_ = 0.0;
    double xup_ = 100.0;
    bool has_binning_ = false;

    // Binning (Y)
    int nbinsy_ = 100;
    double ylow_ = 0.0;
    double yup_ = 100.0;
    bool has_binning_y_ = false;

    // Binning (Z)
    int nbinsz_ = 100;
    double zlow_ = 0.0;
    double zup_ = 100.0;
    bool has_binning_z_ = false;

    // Metadata
    std::string folder_;
    std::string description_;
    std::vector<std::string> tags_;
};

// ============================================================================
// Convenience function for starting builder
// ============================================================================

/**
 * @brief Convenience function to start building a histogram
 *
 * Example:
 *   auto hist = histogram()
 *                   .name("h_p")
 *                   .bins(100, 0, 3000)
 *                   .build1D();
 */
inline HistogramBuilder histogram() {
    return HistogramBuilder();
}

#endif // HISTOGRAM_BUILDER_H
