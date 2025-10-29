/**
 * @file HistogramFactory.h
 * @brief Factory for creating histograms and histogram arrays
 *
 * This file provides factory methods for creating histograms with reduced boilerplate.
 * It replaces repetitive sprintf loops with clean, expressive factory methods.
 *
 * Key features:
 * - Create single histograms with sensible defaults
 * - Create 1D arrays of histograms (e.g., h_p[10])
 * - Create 2D matrices of histograms (e.g., h_theta[10][5])
 * - Automatic naming with indices
 * - Integration with HistogramRegistry
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef HISTOGRAM_FACTORY_H
#define HISTOGRAM_FACTORY_H

#include <string>
#include <vector>
#include <memory>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

// ============================================================================
// HistogramFactory: Factory methods for histogram creation
// ============================================================================

class HistogramFactory {
public:

    // ========================================================================
    // 1D Histograms (TH1F)
    // ========================================================================

    /**
     * @brief Create a single 1D histogram
     *
     * @param name Histogram name
     * @param title Histogram title (use empty string for same as name)
     * @param nbins Number of bins
     * @param xlow Lower edge of first bin
     * @param xup Upper edge of last bin
     * @return std::unique_ptr<TH1F> Ownership of created histogram
     */
    static std::unique_ptr<TH1F> create1D(
        const std::string& name,
        const std::string& title,
        int nbins,
        double xlow,
        double xup)
    {
        const char* title_cstr = title.empty() ? name.c_str() : title.c_str();
        return std::make_unique<TH1F>(name.c_str(), title_cstr, nbins, xlow, xup);
    }

    /**
     * @brief Create a 1D histogram array
     *
     * Creates histograms named: basename_0, basename_1, ..., basename_(n-1)
     *
     * Example:
     *   auto hists = create1DArray("h_theta", "Theta", 10, 100, 0, 180);
     *   // Creates: h_theta_0, h_theta_1, ..., h_theta_9
     *
     * @param basename Base name for histograms
     * @param base_title Base title (will append index)
     * @param array_size Number of histograms to create
     * @param nbins Number of bins per histogram
     * @param xlow Lower edge
     * @param xup Upper edge
     * @return std::vector<std::unique_ptr<TH1F>> Vector of histograms
     */
    static std::vector<std::unique_ptr<TH1F>> create1DArray(
        const std::string& basename,
        const std::string& base_title,
        int array_size,
        int nbins,
        double xlow,
        double xup)
    {
        std::vector<std::unique_ptr<TH1F>> histograms;
        histograms.reserve(array_size);

        for (int i = 0; i < array_size; ++i) {
            std::string name = basename + "_" + std::to_string(i);
            std::string title = base_title + " [" + std::to_string(i) + "]";

            histograms.push_back(
                std::make_unique<TH1F>(name.c_str(), title.c_str(), nbins, xlow, xup)
            );
        }

        return histograms;
    }

    /**
     * @brief Create a 2D matrix of 1D histograms
     *
     * Creates histograms named: basename_i_j
     *
     * Example:
     *   auto matrix = create1DMatrix("h_p", "Momentum", 10, 5, 100, 0, 3000);
     *   // Creates: h_p_0_0, h_p_0_1, ..., h_p_9_4
     *
     * @param basename Base name
     * @param base_title Base title
     * @param rows Number of rows
     * @param cols Number of columns
     * @param nbins Number of bins per histogram
     * @param xlow Lower edge
     * @param xup Upper edge
     * @return 2D vector of histograms
     */
    static std::vector<std::vector<std::unique_ptr<TH1F>>> create1DMatrix(
        const std::string& basename,
        const std::string& base_title,
        int rows,
        int cols,
        int nbins,
        double xlow,
        double xup)
    {
        std::vector<std::vector<std::unique_ptr<TH1F>>> matrix;
        matrix.reserve(rows);

        for (int i = 0; i < rows; ++i) {
            std::vector<std::unique_ptr<TH1F>> row;
            row.reserve(cols);

            for (int j = 0; j < cols; ++j) {
                std::string name = basename + "_" + std::to_string(i) + "_" + std::to_string(j);
                std::string title = base_title + " [" + std::to_string(i) + "][" + std::to_string(j) + "]";

                row.push_back(
                    std::make_unique<TH1F>(name.c_str(), title.c_str(), nbins, xlow, xup)
                );
            }

            matrix.push_back(std::move(row));
        }

        return matrix;
    }

    // ========================================================================
    // 2D Histograms (TH2F)
    // ========================================================================

    /**
     * @brief Create a single 2D histogram
     */
    static std::unique_ptr<TH2F> create2D(
        const std::string& name,
        const std::string& title,
        int nbinsx, double xlow, double xup,
        int nbinsy, double ylow, double yup)
    {
        const char* title_cstr = title.empty() ? name.c_str() : title.c_str();
        return std::make_unique<TH2F>(
            name.c_str(), title_cstr,
            nbinsx, xlow, xup,
            nbinsy, ylow, yup
        );
    }

    /**
     * @brief Create a 2D histogram array
     */
    static std::vector<std::unique_ptr<TH2F>> create2DArray(
        const std::string& basename,
        const std::string& base_title,
        int array_size,
        int nbinsx, double xlow, double xup,
        int nbinsy, double ylow, double yup)
    {
        std::vector<std::unique_ptr<TH2F>> histograms;
        histograms.reserve(array_size);

        for (int i = 0; i < array_size; ++i) {
            std::string name = basename + "_" + std::to_string(i);
            std::string title = base_title + " [" + std::to_string(i) + "]";

            histograms.push_back(
                std::make_unique<TH2F>(
                    name.c_str(), title.c_str(),
                    nbinsx, xlow, xup,
                    nbinsy, ylow, yup
                )
            );
        }

        return histograms;
    }

    /**
     * @brief Create a 2D matrix of 2D histograms
     */
    static std::vector<std::vector<std::unique_ptr<TH2F>>> create2DMatrix(
        const std::string& basename,
        const std::string& base_title,
        int rows,
        int cols,
        int nbinsx, double xlow, double xup,
        int nbinsy, double ylow, double yup)
    {
        std::vector<std::vector<std::unique_ptr<TH2F>>> matrix;
        matrix.reserve(rows);

        for (int i = 0; i < rows; ++i) {
            std::vector<std::unique_ptr<TH2F>> row;
            row.reserve(cols);

            for (int j = 0; j < cols; ++j) {
                std::string name = basename + "_" + std::to_string(i) + "_" + std::to_string(j);
                std::string title = base_title + " [" + std::to_string(i) + "][" + std::to_string(j) + "]";

                row.push_back(
                    std::make_unique<TH2F>(
                        name.c_str(), title.c_str(),
                        nbinsx, xlow, xup,
                        nbinsy, ylow, yup
                    )
                );
            }

            matrix.push_back(std::move(row));
        }

        return matrix;
    }

    // ========================================================================
    // 3D Histograms (TH3F)
    // ========================================================================

    /**
     * @brief Create a single 3D histogram
     */
    static std::unique_ptr<TH3F> create3D(
        const std::string& name,
        const std::string& title,
        int nbinsx, double xlow, double xup,
        int nbinsy, double ylow, double yup,
        int nbinsz, double zlow, double zup)
    {
        const char* title_cstr = title.empty() ? name.c_str() : title.c_str();
        return std::make_unique<TH3F>(
            name.c_str(), title_cstr,
            nbinsx, xlow, xup,
            nbinsy, ylow, yup,
            nbinsz, zlow, zup
        );
    }

    /**
     * @brief Create a 3D histogram array
     */
    static std::vector<std::unique_ptr<TH3F>> create3DArray(
        const std::string& basename,
        const std::string& base_title,
        int array_size,
        int nbinsx, double xlow, double xup,
        int nbinsy, double ylow, double yup,
        int nbinsz, double zlow, double zup)
    {
        std::vector<std::unique_ptr<TH3F>> histograms;
        histograms.reserve(array_size);

        for (int i = 0; i < array_size; ++i) {
            std::string name = basename + "_" + std::to_string(i);
            std::string title = base_title + " [" + std::to_string(i) + "]";

            histograms.push_back(
                std::make_unique<TH3F>(
                    name.c_str(), title.c_str(),
                    nbinsx, xlow, xup,
                    nbinsy, ylow, yup,
                    nbinsz, zlow, zup
                )
            );
        }

        return histograms;
    }

    // ========================================================================
    // Convenience: Create and register in one call
    // ========================================================================

    /**
     * @brief Create 1D histogram and add to registry
     *
     * Example:
     *   HistogramFactory::createAndRegister1D(
     *       registry, "h_theta", "Theta", 100, 0, 180, "angular"
     *   );
     */
    template<typename Registry>
    static void createAndRegister1D(
        Registry& registry,
        const std::string& name,
        const std::string& title,
        int nbins, double xlow, double xup,
        const std::string& folder = "",
        const std::string& description = "")
    {
        auto hist = create1D(name, title, nbins, xlow, xup);
        registry.add(std::move(hist), folder, description);
    }

    /**
     * @brief Create 1D array and add all to registry
     *
     * Example:
     *   HistogramFactory::createAndRegister1DArray(
     *       registry, "h_theta", "Theta", 10, 100, 0, 180, "angular"
     *   );
     */
    template<typename Registry>
    static void createAndRegister1DArray(
        Registry& registry,
        const std::string& basename,
        const std::string& base_title,
        int array_size,
        int nbins, double xlow, double xup,
        const std::string& folder = "",
        const std::string& description = "")
    {
        auto histograms = create1DArray(basename, base_title, array_size, nbins, xlow, xup);

        for (auto& hist : histograms) {
            std::string name = hist->GetName();
            registry.add(std::move(hist), folder, description);
        }
    }

    /**
     * @brief Create 2D histogram and add to registry
     */
    template<typename Registry>
    static void createAndRegister2D(
        Registry& registry,
        const std::string& name,
        const std::string& title,
        int nbinsx, double xlow, double xup,
        int nbinsy, double ylow, double yup,
        const std::string& folder = "",
        const std::string& description = "")
    {
        auto hist = create2D(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);
        registry.add(std::move(hist), folder, description);
    }

    /**
     * @brief Create 2D array and add all to registry
     */
    template<typename Registry>
    static void createAndRegister2DArray(
        Registry& registry,
        const std::string& basename,
        const std::string& base_title,
        int array_size,
        int nbinsx, double xlow, double xup,
        int nbinsy, double ylow, double yup,
        const std::string& folder = "",
        const std::string& description = "")
    {
        auto histograms = create2DArray(basename, base_title, array_size,
                                       nbinsx, xlow, xup,
                                       nbinsy, ylow, yup);

        for (auto& hist : histograms) {
            std::string name = hist->GetName();
            registry.add(std::move(hist), folder, description);
        }
    }

private:
    // Static class - no instances
    HistogramFactory() = delete;
};

#endif // HISTOGRAM_FACTORY_H
