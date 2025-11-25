/**
 * @file HistogramRegistry.h
 * @brief Centralized registry for histogram and ntuple management
 *
 * This file provides a modern histogram management system that replaces global
 * histogram pointers with a centralized registry. Key features:
 *
 * - RAII memory management using std::unique_ptr
 * - Named access to histograms (no more global pointers)
 * - Metadata support (folders, descriptions, tags)
 * - Unified storage for TH1F, TH2F, TH3F via polymorphism
 * - Support for HNtuple objects
 * - Automatic ROOT file organization (folders)
 * - Query capabilities (list by folder, search by tag)
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef HISTOGRAM_REGISTRY_H
#define HISTOGRAM_REGISTRY_H

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TDirectory.h>
#include "hntuple.h"

// ============================================================================
// HistogramMetadata: Store information about each histogram
// ============================================================================

struct HistogramMetadata {
    std::string name;           // Unique identifier
    std::string folder;         // ROOT folder path (e.g., "proton/cms")
    std::string description;    // Human-readable description
    std::vector<std::string> tags;  // Tags for searching (e.g., "angular", "momentum")

    HistogramMetadata() = default;

    HistogramMetadata(const std::string& n,
                     const std::string& f = "",
                     const std::string& d = "")
        : name(n), folder(f), description(d) {}

    // Add tag
    void addTag(const std::string& tag) {
        tags.push_back(tag);
    }

    // Check if has tag
    bool hasTag(const std::string& tag) const {
        return std::find(tags.begin(), tags.end(), tag) != tags.end();
    }
};

// ============================================================================
// HistogramRegistry: Centralized histogram storage and management
// ============================================================================

class HistogramRegistry {
public:

    // ------------------------------------------------------------------------
    // Singleton access (optional - can also be used as regular object)
    // ------------------------------------------------------------------------
    static HistogramRegistry& instance() {
        static HistogramRegistry registry;
        return registry;
    }

    // ------------------------------------------------------------------------
    // Add histogram to registry (takes ownership via unique_ptr)
    // ------------------------------------------------------------------------
    void add(std::unique_ptr<TH1> hist, const HistogramMetadata& meta) {
        if (!hist) {
            throw std::runtime_error("HistogramRegistry::add() - Cannot add null histogram!");
        }

        const std::string& name = meta.name;

        // Check for duplicates
        if (histograms_.find(name) != histograms_.end()) {
            throw std::runtime_error("HistogramRegistry::add() - Histogram '" +
                                   name + "' already exists!");
        }

        // Store histogram and metadata
        histograms_[name] = std::move(hist);
        metadata_[name] = meta;
    }

    // ------------------------------------------------------------------------
    // Add histogram with automatic metadata from histogram name/title
    // ------------------------------------------------------------------------
    void add(std::unique_ptr<TH1> hist,
             const std::string& folder = "",
             const std::string& description = "") {
        if (!hist) {
            throw std::runtime_error("HistogramRegistry::add() - Cannot add null histogram!");
        }

        HistogramMetadata meta;
        meta.name = hist->GetName();
        meta.folder = folder;
        meta.description = description.empty() ? hist->GetTitle() : description;

        add(std::move(hist), meta);
    }

    // ------------------------------------------------------------------------
    // Add HNtuple to registry
    // ------------------------------------------------------------------------
    void addNtuple(std::unique_ptr<HNtuple> ntuple,
                   const std::string& folder = "",
                   const std::string& description = "") {
        if (!ntuple) {
            throw std::runtime_error("HistogramRegistry::addNtuple() - Cannot add null ntuple!");
        }

        // Copy strings before moving the ntuple
        std::string name(ntuple->getName());
        std::string title(ntuple->getTitle());

        if (ntuples_.find(name) != ntuples_.end()) {
            throw std::runtime_error("HistogramRegistry::addNtuple() - Ntuple '" +
                                   name + "' already exists!");
        }

        NtupleMetadata meta;
        meta.name = name;
        meta.folder = folder;
        meta.description = description.empty() ? title : description;

        ntuples_[name] = std::move(ntuple);
        ntuple_metadata_[name] = meta;
    }

    // ------------------------------------------------------------------------
    // Get histogram (non-const version for filling)
    // ------------------------------------------------------------------------
    TH1* get(const std::string& name) {
        auto it = histograms_.find(name);
        if (it == histograms_.end()) {
            throw std::runtime_error("HistogramRegistry::get() - Histogram '" +
                                   name + "' not found!");
        }
        return it->second.get();
    }

    // ------------------------------------------------------------------------
    // Get histogram (const version for reading)
    // ------------------------------------------------------------------------
    const TH1* get(const std::string& name) const {
        auto it = histograms_.find(name);
        if (it == histograms_.end()) {
            throw std::runtime_error("HistogramRegistry::get() - Histogram '" +
                                   name + "' not found!");
        }
        return it->second.get();
    }

    // ------------------------------------------------------------------------
    // Get histogram with type checking (more convenient)
    // ------------------------------------------------------------------------
    template<typename T>
    T* getAs(const std::string& name) {
        TH1* hist = get(name);
        T* typed = dynamic_cast<T*>(hist);
        if (!typed) {
            throw std::runtime_error("HistogramRegistry::getAs() - Histogram '" +
                                   name + "' is not of requested type!");
        }
        return typed;
    }

    // ------------------------------------------------------------------------
    // Get HNtuple
    // ------------------------------------------------------------------------
    HNtuple* getNtuple(const std::string& name) {
        auto it = ntuples_.find(name);
        if (it == ntuples_.end()) {
            throw std::runtime_error("HistogramRegistry::getNtuple() - Ntuple '" +
                                   name + "' not found!");
        }
        return it->second.get();
    }

    // ------------------------------------------------------------------------
    // Check if histogram exists
    // ------------------------------------------------------------------------
    bool has(const std::string& name) const {
        return histograms_.find(name) != histograms_.end();
    }

    // ------------------------------------------------------------------------
    // Check if ntuple exists
    // ------------------------------------------------------------------------
    bool hasNtuple(const std::string& name) const {
        return ntuples_.find(name) != ntuples_.end();
    }

    // ------------------------------------------------------------------------
    // List all histogram names
    // ------------------------------------------------------------------------
    std::vector<std::string> listAll() const {
        std::vector<std::string> names;
        names.reserve(histograms_.size());
        for (const auto& pair : histograms_) {
            names.push_back(pair.first);
        }
        return names;
    }

    // ------------------------------------------------------------------------
    // List histograms in a specific folder
    // ------------------------------------------------------------------------
    std::vector<std::string> listByFolder(const std::string& folder) const {
        std::vector<std::string> names;
        for (const auto& pair : metadata_) {
            if (pair.second.folder == folder) {
                names.push_back(pair.first);
            }
        }
        return names;
    }

    // ------------------------------------------------------------------------
    // List histograms with a specific tag
    // ------------------------------------------------------------------------
    std::vector<std::string> listByTag(const std::string& tag) const {
        std::vector<std::string> names;
        for (const auto& pair : metadata_) {
            if (pair.second.hasTag(tag)) {
                names.push_back(pair.first);
            }
        }
        return names;
    }

    // ------------------------------------------------------------------------
    // Get metadata for a histogram
    // ------------------------------------------------------------------------
    const HistogramMetadata& getMetadata(const std::string& name) const {
        auto it = metadata_.find(name);
        if (it == metadata_.end()) {
            throw std::runtime_error("HistogramRegistry::getMetadata() - Histogram '" +
                                   name + "' not found!");
        }
        return it->second;
    }

    // ------------------------------------------------------------------------
    // Write all histograms to file (organized by folders)
    // ------------------------------------------------------------------------
    void writeToFile(TFile* file) {
        if (!file || !file->IsOpen()) {
            throw std::runtime_error("HistogramRegistry::writeToFile() - File is not open!");
        }

        // Group histograms by folder
        std::map<std::string, std::vector<std::string>> folder_contents;

        for (const auto& pair : metadata_) {
            const std::string& name = pair.first;
            const std::string& folder = pair.second.folder;

            if (folder.empty()) {
                // Write to root directory
                file->cd();
                histograms_[name]->Write();
                // Release ownership - ROOT now owns this histogram
                histograms_[name].release();
            } else {
                folder_contents[folder].push_back(name);
            }
        }

        // Create folders and write histograms
        for (const auto& folder_pair : folder_contents) {
            const std::string& folder = folder_pair.first;
            const auto& hist_names = folder_pair.second;

            // Create folder hierarchy (e.g., "proton/cms" creates "proton" then "cms")
            TDirectory* dir = createFolderHierarchy(file, folder);
            dir->cd();

            // Write all histograms in this folder
            for (const auto& name : hist_names) {
                histograms_[name]->Write();
                // Release ownership - ROOT now owns this histogram
                histograms_[name].release();
            }
        }

        // Write ntuples
        for (auto& pair : ntuples_) {
            const std::string& name = pair.first;
            const auto& meta = ntuple_metadata_[name];

            if (meta.folder.empty()) {
                file->cd();
            } else {
                TDirectory* dir = createFolderHierarchy(file, meta.folder);
                dir->cd();
            }

            pair.second->Write();
            // Release ownership - ROOT now owns this ntuple
            pair.second.release();
        }

        // Return to root directory
        file->cd();
    }

    // ------------------------------------------------------------------------
    // Print summary of registry contents
    // ------------------------------------------------------------------------
    void printSummary(std::ostream& os = std::cout) const {
        os << "╔════════════════════════════════════════════════════════════════╗\n";
        os << "║  HistogramRegistry Summary                                     ║\n";
        os << "╠════════════════════════════════════════════════════════════════╣\n";
        // Dynamic padding for Total histograms
        os << "║ Total histograms: " << histograms_.size() 
           << std::string(45 - std::to_string(histograms_.size()).length(), ' ') << "║\n";
        // Dynamic padding for Total ntuples
        os << "║ Total ntuples:    " << ntuples_.size() 
           << std::string(45 - std::to_string(ntuples_.size()).length(), ' ') << "║\n";

        // Count by folder
        std::map<std::string, int> folder_counts;
        for (const auto& pair : metadata_) {
            const std::string& folder = pair.second.folder;
            folder_counts[folder.empty() ? "[root]" : folder]++;
        }

        os << "║                                                                ║\n";
        os << "║ Histograms by folder:                                          ║\n";
        for (const auto& fc : folder_counts) {
            // Calculate exact padding: 66 total chars - 1(║) - 3(spaces) - folder - 2(": ") - number - 1(║)
            // = 66 - 7 - folder - number = 59 - folder - number
            int padding = 59 - fc.first.length() - std::to_string(fc.second).length();
            os << "║   " << fc.first << ": " << fc.second 
               << std::string(padding, ' ') << "║\n";
        }
        os << "╚════════════════════════════════════════════════════════════════╝\n";
    }

    // ------------------------------------------------------------------------
    // Get statistics
    // ------------------------------------------------------------------------
    size_t size() const { return histograms_.size(); }
    size_t ntupleCount() const { return ntuples_.size(); }
    bool empty() const { return histograms_.empty() && ntuples_.empty(); }

    // ------------------------------------------------------------------------
    // Clear all histograms (useful for testing)
    // ------------------------------------------------------------------------
    void clear() {
        histograms_.clear();
        metadata_.clear();
        ntuples_.clear();
        ntuple_metadata_.clear();
    }

private:
    // Metadata for ntuples
    struct NtupleMetadata {
        std::string name;
        std::string folder;
        std::string description;
    };

    // Storage
    std::map<std::string, std::unique_ptr<TH1>> histograms_;
    std::map<std::string, HistogramMetadata> metadata_;
    std::map<std::string, std::unique_ptr<HNtuple>> ntuples_;
    std::map<std::string, NtupleMetadata> ntuple_metadata_;

    // Helper: Create folder hierarchy in ROOT file
    TDirectory* createFolderHierarchy(TFile* file, const std::string& path) const {
        TDirectory* current = file;

        // Split path by '/'
        std::vector<std::string> parts;
        std::string part;
        for (char c : path) {
            if (c == '/') {
                if (!part.empty()) {
                    parts.push_back(part);
                    part.clear();
                }
            } else {
                part += c;
            }
        }
        if (!part.empty()) {
            parts.push_back(part);
        }

        // Create each folder
        for (const auto& folder_name : parts) {
            TDirectory* subdir = current->GetDirectory(folder_name.c_str());
            if (!subdir) {
                subdir = current->mkdir(folder_name.c_str());
            }
            current = subdir;
        }

        return current;
    }

public:
    // Default constructor - public for non-singleton usage
    HistogramRegistry() = default;
    
    // Delete copy/move for safety
    HistogramRegistry(const HistogramRegistry&) = delete;
    HistogramRegistry& operator=(const HistogramRegistry&) = delete;
};

#endif // HISTOGRAM_REGISTRY_H
