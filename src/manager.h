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
#include "dynamic_hntuple.h"
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
        // Finalize any dynamic ntuples before closing
        for (auto& pair : dynamic_ntuples_) {
            if (!pair.second->isFinalized()) {
                pair.second->finalize();
            }
        }
        
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

        // Finalize all dynamic ntuples (TTree → TNtuple conversion)
        for (auto& pair : dynamic_ntuples_) {
            if (!pair.second->isFinalized()) {
                pair.second->finalize();
            }
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
     * @brief Create and register ntuple (original HNtuple - requires prebooking)
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
    // DynamicHNtuple management (RECOMMENDED - no prebooking needed!)
    // ------------------------------------------------------------------------

    /**
     * @brief Create DynamicHNtuple - add variables at ANY time!
     *
     * Unlike regular HNtuple, DynamicHNtuple allows adding variables at any
     * point during processing. The structure is finalized at closeFile().
     *
     * Example:
     *   auto& nt = manager.createDynamicNtuple("events", "Event data");
     *   // In event loop - add ANY variable at ANY time:
     *   nt["mass"] = mass_value;
     *   nt["event"] = event_num;
     *   if (has_rare_condition) {
     *       nt["rare_var"] = value;  // OK even if first time!
     *   }
     *   nt.fill();
     *
     * @param name Name of the ntuple
     * @param title Title/description
     * @param missing_value Value for variables not set in an event (default: -1)
     * @param keep_intermediate Keep intermediate TTree file (default: false)
     * @return Reference to DynamicHNtuple for direct access
     */
    DynamicHNtuple& createDynamicNtuple(const std::string& name,
                                        const std::string& title = "",
                                        Float_t missing_value = -1.0f,
                                        bool keep_intermediate = false)
    {
        if (!file_ || !file_->IsOpen()) {
            throw std::runtime_error("Manager::createDynamicNtuple() - No file open! Call openFile() first.");
        }

        if (dynamic_ntuples_.find(name) != dynamic_ntuples_.end()) {
            throw std::runtime_error("Manager::createDynamicNtuple() - DynamicNtuple '" + 
                                   name + "' already exists!");
        }

        auto ntuple = std::make_unique<DynamicHNtuple>(
            name, 
            title.empty() ? name : title,
            file_.get(),
            missing_value,
            keep_intermediate
        );

        dynamic_ntuples_[name] = std::move(ntuple);
        return *dynamic_ntuples_[name];
    }

    /**
     * @brief Get DynamicHNtuple by name (reference for clean syntax)
     *
     * Example:
     *   auto& nt = manager.getDynamicNtuple("events");
     *   nt["var"] = value;  // Clean syntax without (*nt)["var"]
     */
    DynamicHNtuple& getDynamicNtuple(const std::string& name) {
        auto it = dynamic_ntuples_.find(name);
        if (it == dynamic_ntuples_.end()) {
            throw std::runtime_error("Manager::getDynamicNtuple() - DynamicNtuple '" +
                                   name + "' not found!");
        }
        return *(it->second);
    }

    /**
     * @brief Check if DynamicHNtuple exists
     */
    bool hasDynamicNtuple(const std::string& name) const {
        return dynamic_ntuples_.find(name) != dynamic_ntuples_.end();
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
    // PT3 / PT2 trigger flag — drives auto-routing in fill().
    //
    // Each event in main.cc is classified by trigbit:
    //     trigbit == 8192  →  TRIG_PT3   (biased dilepton trigger; main analysis)
    //     trigbit == 4096  →  TRIG_PT2   (downscaled unbiased trigger;
    //                                     control sample for trigger correction)
    //     anything else    →  TRIG_NONE  (event skipped)
    //
    // setEventTriggerType() is called ONCE per event before processEvent.
    // fill()/fillw() then route to either H (PT3) or H_pt2 (PT2). Both
    // versions of every histogram are created by createPT2Clones() right
    // after setupHistograms — so downstream analysis can compute the
    // per-bin trigger correction directly as PT2/PT3 ratio (no time-segment
    // averaging buried inside the value).
    // ------------------------------------------------------------------------
    enum TriggerType { TRIG_NONE, TRIG_PT3, TRIG_PT2 };

    /// Set the per-event trigger flag. In addition to driving histogram
    /// auto-routing in fill(), this stamps a `trigbit` branch on every
    /// dynamic ntuple so downstream analyses can filter by trigger:
    ///     trigbit == 8192  → PT3
    ///     trigbit == 4096  → PT2
    /// The stamp happens here (once per event) and persists through any
    /// subsequent nt.fill() because DynamicHNtuple only resets fields it
    /// has written — `trigbit` stays committed for the row.
    void setEventTriggerType(TriggerType t) {
        current_trigger_ = t;
        const Float_t code = (t == TRIG_PT3) ? 8192.0f
                           : (t == TRIG_PT2) ? 4096.0f
                           : 0.0f;
        for (auto& pair : dynamic_ntuples_) {
            if (pair.second && !pair.second->isFinalized()) {
                (*pair.second)["trigbit"] = code;
            }
        }
    }
    TriggerType eventTriggerType() const { return current_trigger_; }

    /// Clone every registered histogram with a "_pt2" suffix (same folder,
    /// same metadata, empty contents). Call once after setupHistograms.
    /// Idempotent: existing "_pt2" entries are left alone.
    void createPT2Clones() {
        std::vector<std::string> names = registry_.listAll();
        for (const auto& name : names) {
            const std::string suffix = "_pt2";
            if (name.size() >= suffix.size() &&
                name.compare(name.size() - suffix.size(),
                             suffix.size(), suffix) == 0) {
                continue;  // already a _pt2 clone
            }
            const std::string pt2_name = name + suffix;
            if (registry_.has(pt2_name)) continue;

            TH1* src = registry_.get(name);
            std::unique_ptr<TH1> clone(
                static_cast<TH1*>(src->Clone(pt2_name.c_str())));
            clone->Reset();
            // Append " (PT2)" to the title portion (the part before the
            // first ';' which separates title from axis labels in ROOT).
            std::string t = src->GetTitle();
            const size_t semi = t.find(';');
            if (semi == std::string::npos) {
                t += " (PT2)";
            } else {
                t.insert(semi, " (PT2)");
            }
            clone->SetTitle(t.c_str());
            clone->SetDirectory(nullptr);  // registry owns; not gDirectory

            HistogramMetadata meta = registry_.getMetadata(name);
            meta.name = pt2_name;
            registry_.add(std::move(clone), meta);
        }
    }

    // ------------------------------------------------------------------------
    // Fill helpers (PT3/PT2 auto-routed).
    //
    //   current_trigger_ == TRIG_PT3 → fill <name>
    //   current_trigger_ == TRIG_PT2 → fill <name>_pt2
    //   current_trigger_ == TRIG_NONE → no-op (defensive)
    // ------------------------------------------------------------------------
private:
    const std::string& targetName_(const std::string& name) const {
        if (current_trigger_ == TRIG_PT2) {
            // Resolve _pt2 target name lazily; static thread-unsafe storage
            // is fine because main.cc is single-threaded.
            static thread_local std::string buf;
            buf = name + "_pt2";
            return buf;
        }
        return name;
    }

public:
    /// Fill 1D histogram (unweighted), auto-routed by current trigger flag.
    void fill(const std::string& name, double value) {
        if (current_trigger_ == TRIG_NONE) return;
        getHistogram(targetName_(name))->Fill(value);
    }

    /// Fill 1D histogram with explicit weight (auto-routed).
    void fillw(const std::string& name, double value, double weight) {
        if (current_trigger_ == TRIG_NONE) return;
        getHistogram(targetName_(name))->Fill(value, weight);
    }

    /// Fill 2D histogram (unweighted), auto-routed.
    void fill(const std::string& name, double x, double y) {
        if (current_trigger_ == TRIG_NONE) return;
        getHistogramAs<TH2>(targetName_(name))->Fill(x, y);
    }

    /// Fill 2D histogram with explicit weight, auto-routed.
    void fillw(const std::string& name, double x, double y, double weight) {
        if (current_trigger_ == TRIG_NONE) return;
        getHistogramAs<TH2>(targetName_(name))->Fill(x, y, weight);
    }

    /// Fill 3D histogram (unweighted), auto-routed.
    void fill(const std::string& name, double x, double y, double z) {
        if (current_trigger_ == TRIG_NONE) return;
        getHistogramAs<TH3>(targetName_(name))->Fill(x, y, z);
    }

    /// Fill 3D histogram with explicit weight, auto-routed.
    void fillw(const std::string& name, double x, double y, double z, double weight) {
        if (current_trigger_ == TRIG_NONE) return;
        getHistogramAs<TH3>(targetName_(name))->Fill(x, y, z, weight);
    }

    // ------------------------------------------------------------------------
    // Statistics and diagnostics
    // ------------------------------------------------------------------------

    /**
     * @brief Print summary of all managed histograms/ntuples
     */
    void printSummary(std::ostream& os = std::cout) const {
        registry_.printSummary(os, dynamic_ntuples_.size());
    }

    /**
     * @brief Get number of histograms
     */
    size_t histogramCount() const {
        return registry_.size();
    }

    /**
     * @brief Get number of ntuples (HNtuple only)
     */
    size_t ntupleCount() const {
        return registry_.ntupleCount();
    }

    /**
     * @brief Get number of dynamic ntuples
     */
    size_t dynamicNtupleCount() const {
        return dynamic_ntuples_.size();
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

    // Dynamic ntuples (managed separately due to finalization needs)
    std::map<std::string, std::unique_ptr<DynamicHNtuple>> dynamic_ntuples_;

    // Current event's trigger classification — set once per event by main.cc
    // before processEvent. Drives auto-routing of fill()/fillw() to either
    // the H (PT3) or the H_pt2 (PT2) variant of every registered histogram.
    TriggerType current_trigger_ = TRIG_NONE;
};

#endif // MANAGER_H
