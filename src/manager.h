/**
 * @file Manager.h
 * @brief Modern histogram and ntuple manager using modern C++ patterns
 *
 * Key features:
 * - Uses HistogramRegistry for centralized storage
 * - Uses smart pointers (RAII) instead of raw pointers
 * - Eliminates global histogram pointers
 * - Provides type-safe access
 * - Supports metadata and folder organization
 * - Fixes memory leaks from original implementation
 *
 * This class is designed to be backward-compatible with existing code while
 * providing a migration path to the new architecture.
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef MANAGER_H
#define MANAGER_H

#include <string>
#include <memory>
#include <stdexcept>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include "hntuple.h"
#include "histogram_registry.h"
#include "histogram_factory.h"
#include "histogram_builder.h"

// ============================================================================
// Manager: Modern histogram/ntuple manager
// ============================================================================

class Manager {
public:

    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    Manager() = default;

    // ------------------------------------------------------------------------
    // Destructor (RAII - automatic cleanup)
    // ------------------------------------------------------------------------
    ~Manager() {
        // Registry uses smart pointers, so cleanup is automatic!
        // No memory leaks (unlike original Manager with commented-out delete)
        if (file_ && file_->IsOpen()) {
            file_->Close();
        }
    }

    // Delete copy and move operations (registry is not movable)
    Manager(const Manager&) = delete;
    Manager& operator=(const Manager&) = delete;
    Manager(Manager&&) = delete;
    Manager& operator=(Manager&&) = delete;

    // ------------------------------------------------------------------------
    // File management
    // ------------------------------------------------------------------------

    /**
     * @brief Open output ROOT file
     *
     * @param filename Path to output file
     * @param option ROOT file option ("RECREATE", "UPDATE", etc.)
     */
    void openFile(const std::string& filename, const std::string& option = "RECREATE") {
        if (file_ && file_->IsOpen()) {
            throw std::runtime_error("Manager::openFile() - File already open!");
        }

        file_ = std::make_unique<TFile>(filename.c_str(), option.c_str());

        if (!file_->IsOpen()) {
            throw std::runtime_error("Manager::openFile() - Failed to open file: " + filename);
        }

        std::cout << "Manager: Opened file '" << filename << "' with option '" << option << "'\n";
    }

    /**
     * @brief Close output file (writes all histograms automatically)
     */
    void closeFile() {
        if (!file_ || !file_->IsOpen()) {
            throw std::runtime_error("Manager::closeFile() - No file is open!");
        }

        // Write all histograms to file with folder organization
        registry_.writeToFile(file_.get());

        file_->Close();
        std::cout << "Manager: File closed and all histograms written.\n";
    }

    /**
     * @brief Get pointer to ROOT file (for direct access if needed)
     */
    TFile* getFile() {
        return file_.get();
    }

    // ------------------------------------------------------------------------
    // Access to registry
    // ------------------------------------------------------------------------

    /**
     * @brief Get reference to histogram registry
     */
    HistogramRegistry& registry() {
        return registry_;
    }

    const HistogramRegistry& registry() const {
        return registry_;
    }

    // ------------------------------------------------------------------------
    // Histogram creation (convenience wrappers)
    // ------------------------------------------------------------------------

    /**
     * @brief Create and register 1D histogram
     *
     * Example:
     *   manager.create1D("h_theta", "Theta", 100, 0, 180, "angular");
     */
    void create1D(const std::string& name,
                  const std::string& title,
                  int nbins, double xlow, double xup,
                  const std::string& folder = "",
                  const std::string& description = "")
    {
        HistogramFactory::createAndRegister1D(
            registry_, name, title, nbins, xlow, xup, folder, description
        );
    }

    /**
     * @brief Create and register 1D histogram array
     *
     * Example:
     *   manager.create1DArray("h_theta", "Theta", 10, 100, 0, 180, "angular");
     *   // Creates: h_theta_0, h_theta_1, ..., h_theta_9
     */
    void create1DArray(const std::string& basename,
                       const std::string& base_title,
                       int array_size,
                       int nbins, double xlow, double xup,
                       const std::string& folder = "",
                       const std::string& description = "")
    {
        HistogramFactory::createAndRegister1DArray(
            registry_, basename, base_title, array_size, nbins, xlow, xup, folder, description
        );
    }

    /**
     * @brief Create and register 2D histogram
     */
    void create2D(const std::string& name,
                  const std::string& title,
                  int nbinsx, double xlow, double xup,
                  int nbinsy, double ylow, double yup,
                  const std::string& folder = "",
                  const std::string& description = "")
    {
        HistogramFactory::createAndRegister2D(
            registry_, name, title, nbinsx, xlow, xup, nbinsy, ylow, yup, folder, description
        );
    }

    /**
     * @brief Create and register 2D histogram array
     */
    void create2DArray(const std::string& basename,
                       const std::string& base_title,
                       int array_size,
                       int nbinsx, double xlow, double xup,
                       int nbinsy, double ylow, double yup,
                       const std::string& folder = "",
                       const std::string& description = "")
    {
        HistogramFactory::createAndRegister2DArray(
            registry_, basename, base_title, array_size,
            nbinsx, xlow, xup, nbinsy, ylow, yup, folder, description
        );
    }

    // ------------------------------------------------------------------------
    // NTuple management
    // ------------------------------------------------------------------------

    /**
     * @brief Create and register ntuple
     *
     * Example:
     *   manager.createNtuple("events", "Event data", "ntuples");
     */
    void createNtuple(const std::string& name,
                     const std::string& title = "",
                     const std::string& folder = "",
                     Int_t bufsize = 32000)
    {
        if (!file_ || !file_->IsOpen()) {
            throw std::runtime_error("Manager::createNtuple() - No file open! Call openFile() first.");
        }

        auto ntuple = std::make_unique<HNtuple>(name.c_str(), title.c_str(), bufsize);
        ntuple->setFile(file_.get());

        registry_.addNtuple(std::move(ntuple), folder, title);
    }

    // ------------------------------------------------------------------------
    // Histogram/Ntuple access
    // ------------------------------------------------------------------------

    /**
     * @brief Get histogram by name
     *
     * Example:
     *   TH1* h = manager.getHistogram("h_theta");
     *   h->Fill(45.0);
     */
    TH1* getHistogram(const std::string& name) {
        return registry_.get(name);
    }

    /**
     * @brief Get histogram with type checking
     *
     * Example:
     *   TH1F* h = manager.getHistogramAs<TH1F>("h_theta");
     */
    template<typename T>
    T* getHistogramAs(const std::string& name) {
        return registry_.getAs<T>(name);
    }

    /**
     * @brief Get ntuple by name
     *
     * Example:
     *   HNtuple* nt = manager.getNtuple("events");
     *   (*nt)["theta"] = 45.0;
     *   nt->fill();
     */
    HNtuple* getNtuple(const std::string& name) {
        return registry_.getNtuple(name);
    }

    // ------------------------------------------------------------------------
    // Fill helpers (shorthand for common operations)
    // ------------------------------------------------------------------------

    /**
     * @brief Fill 1D histogram (shorthand)
     *
     * Example:
     *   manager.fill("h_theta", 45.0);
     */
    void fill(const std::string& name, double value) {
        getHistogram(name)->Fill(value);
    }

    /**
     * @brief Fill 2D histogram (shorthand)
     */
    void fill(const std::string& name, double x, double y) {
        getHistogramAs<TH2>(name)->Fill(x, y);
    }

    /**
     * @brief Fill 3D histogram (shorthand)
     */
    void fill(const std::string& name, double x, double y, double z) {
        getHistogramAs<TH3>(name)->Fill(x, y, z);
    }

    // ------------------------------------------------------------------------
    // Statistics and diagnostics
    // ------------------------------------------------------------------------

    /**
     * @brief Print summary of all managed histograms/ntuples
     */
    void printSummary(std::ostream& os = std::cout) const {
        registry_.printSummary(os);
    }

    /**
     * @brief Get number of histograms
     */
    size_t histogramCount() const {
        return registry_.size();
    }

    /**
     * @brief Get number of ntuples
     */
    size_t ntupleCount() const {
        return registry_.ntupleCount();
    }

    /**
     * @brief Check if histogram exists
     */
    bool hasHistogram(const std::string& name) const {
        return registry_.has(name);
    }

    /**
     * @brief Check if ntuple exists
     */
    bool hasNtuple(const std::string& name) const {
        return registry_.hasNtuple(name);
    }

    // ------------------------------------------------------------------------
    // Advanced: Batch operations
    // ------------------------------------------------------------------------

    /**
     * @brief List all histogram names in a folder
     */
    std::vector<std::string> listHistogramsInFolder(const std::string& folder) const {
        return registry_.listByFolder(folder);
    }

    /**
     * @brief List all histograms with a specific tag
     */
    std::vector<std::string> listHistogramsByTag(const std::string& tag) const {
        return registry_.listByTag(tag);
    }

private:
    // ROOT file for output
    std::unique_ptr<TFile> file_;

    // Centralized histogram/ntuple storage
    HistogramRegistry registry_;
};

#endif // MANAGER_H
