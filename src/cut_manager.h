/**
 * @file cut_manager.h
 * @brief Unified cut management system for physics analysis
 *
 * Provides a centralized system for managing all types of cuts:
 * - 1D cuts (mass windows, momentum ranges, etc.)
 * - 2D graphical cuts (TCutG)
 * - Trigger/flag selections
 * - Event quality cuts
 *
 * Supports:
 * - Named cut access
 * - Cut statistics tracking
 * - JSON configuration (via AnalysisConfig)
 * - Cut flow analysis
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef CUT_MANAGER_H
#define CUT_MANAGER_H

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <TFile.h>
#include <TCutG.h>
#include <TKey.h>

// ============================================================================
// Cut Types
// ============================================================================

/**
 * @struct RangeCut
 * @brief 1D range cut (min <= x <= max)
 */
struct RangeCut {
    std::string name;
    std::string description;
    double min;
    double max;
    bool active = true;
    
    // Statistics
    mutable Long64_t tested = 0;
    mutable Long64_t passed = 0;
    
    RangeCut() : min(0), max(0) {}
    RangeCut(const std::string& n, double lo, double hi, const std::string& desc = "")
        : name(n), description(desc), min(lo), max(hi) {}
    
    bool pass(double value) const {
        ++tested;
        bool result = (value >= min && value <= max);
        if (result) ++passed;
        return result || !active;
    }
    
    double efficiency() const {
        return tested > 0 ? static_cast<double>(passed) / tested : 0.0;
    }
    
    void reset() { tested = passed = 0; }
};

/**
 * @struct TriggerCut
 * @brief Bit mask trigger selection
 */
struct TriggerCut {
    std::string name;
    std::string description;
    int mask;
    bool require_all = false;  // AND vs OR logic
    bool active = true;
    
    mutable Long64_t tested = 0;
    mutable Long64_t passed = 0;
    
    TriggerCut() : mask(0) {}
    TriggerCut(const std::string& n, int m, bool all = false, const std::string& desc = "")
        : name(n), description(desc), mask(m), require_all(all) {}
    
    bool pass(int trigger) const {
        ++tested;
        bool result;
        if (require_all) {
            result = (trigger & mask) == mask;
        } else {
            result = (trigger & mask) != 0;
        }
        if (result) ++passed;
        return result || !active;
    }
    
    double efficiency() const {
        return tested > 0 ? static_cast<double>(passed) / tested : 0.0;
    }
    
    void reset() { tested = passed = 0; }
};

/**
 * @struct GraphicalCut
 * @brief 2D graphical cut wrapper (TCutG)
 */
struct GraphicalCut {
    std::string name;
    std::string description;
    std::unique_ptr<TCutG> cut;
    bool active = true;
    
    mutable Long64_t tested = 0;
    mutable Long64_t passed = 0;
    
    GraphicalCut() = default;
    GraphicalCut(const std::string& n, TCutG* c, const std::string& desc = "")
        : name(n), description(desc), cut(c) {}
    
    // Move constructor
    GraphicalCut(GraphicalCut&& other) noexcept
        : name(std::move(other.name)), description(std::move(other.description)),
          cut(std::move(other.cut)), active(other.active),
          tested(other.tested), passed(other.passed) {}
    
    // Move assignment
    GraphicalCut& operator=(GraphicalCut&& other) noexcept {
        name = std::move(other.name);
        description = std::move(other.description);
        cut = std::move(other.cut);
        active = other.active;
        tested = other.tested;
        passed = other.passed;
        return *this;
    }
    
    bool pass(double x, double y) const {
        if (!cut) return true;
        ++tested;
        bool result = cut->IsInside(x, y);
        if (result) ++passed;
        return result || !active;
    }
    
    double efficiency() const {
        return tested > 0 ? static_cast<double>(passed) / tested : 0.0;
    }
    
    void reset() { tested = passed = 0; }
};

// ============================================================================
// CutManager: Central cut management
// ============================================================================

class CutManager {
public:
    // ========================================================================
    // Constructor
    // ========================================================================
    
    CutManager() = default;
    
    // ========================================================================
    // 1D Range Cuts
    // ========================================================================
    
    /**
     * @brief Define a 1D range cut
     */
    void defineRangeCut(const std::string& name, double min, double max,
                       const std::string& description = "") {
        if (range_cuts_.find(name) != range_cuts_.end()) {
            std::cerr << "Warning: Overwriting existing cut '" << name << "'\n";
        }
        range_cuts_[name] = RangeCut(name, min, max, description);
    }
    
    /**
     * @brief Test value against named range cut
     */
    bool passRangeCut(const std::string& name, double value) const {
        auto it = range_cuts_.find(name);
        if (it == range_cuts_.end()) {
            throw std::runtime_error("CutManager::passRangeCut() - Cut '" + name + "' not defined!");
        }
        return it->second.pass(value);
    }
    
    /**
     * @brief Get reference to range cut (for modification)
     */
    RangeCut& getRangeCut(const std::string& name) {
        auto it = range_cuts_.find(name);
        if (it == range_cuts_.end()) {
            throw std::runtime_error("CutManager::getRangeCut() - Cut '" + name + "' not defined!");
        }
        return it->second;
    }
    
    // ========================================================================
    // Trigger Cuts
    // ========================================================================
    
    /**
     * @brief Define a trigger cut
     * @param mask Bit mask for trigger selection
     * @param require_all If true, require ALL bits set; if false, require ANY
     */
    void defineTriggerCut(const std::string& name, int mask, bool require_all = false,
                         const std::string& description = "") {
        trigger_cuts_[name] = TriggerCut(name, mask, require_all, description);
    }
    
    /**
     * @brief Test trigger value
     */
    bool passTriggerCut(const std::string& name, int trigger) const {
        auto it = trigger_cuts_.find(name);
        if (it == trigger_cuts_.end()) {
            throw std::runtime_error("CutManager::passTriggerCut() - Cut '" + name + "' not defined!");
        }
        return it->second.pass(trigger);
    }
    
    // ========================================================================
    // 2D Graphical Cuts (TCutG)
    // ========================================================================
    
    /**
     * @brief Load graphical cut from ROOT file
     * @param name Name to register cut under
     * @param filename ROOT file containing TCutG
     * @param cutname Name of TCutG in file (if different from 'name')
     */
    void loadGraphicalCut(const std::string& name, const std::string& filename,
                         const std::string& cutname = "", const std::string& description = "") {
        TFile* f = TFile::Open(filename.c_str(), "READ");
        if (!f || f->IsZombie()) {
            throw std::runtime_error("CutManager::loadGraphicalCut() - Cannot open file: " + filename);
        }
        
        std::string objname = cutname.empty() ? name : cutname;
        TCutG* cut = dynamic_cast<TCutG*>(f->Get(objname.c_str()));
        if (!cut) {
            f->Close();
            throw std::runtime_error("CutManager::loadGraphicalCut() - TCutG '" + objname + 
                                   "' not found in " + filename);
        }
        
        // Clone to own memory (file will close)
        TCutG* cut_clone = dynamic_cast<TCutG*>(cut->Clone());
        cut_clone->SetName(name.c_str());
        f->Close();
        
        graphical_cuts_[name] = GraphicalCut(name, cut_clone, description);
        std::cout << "CutManager: Loaded graphical cut '" << name << "' from " << filename << "\n";
    }
    
    /**
     * @brief Add graphical cut directly
     */
    void addGraphicalCut(const std::string& name, TCutG* cut,
                        const std::string& description = "") {
        if (!cut) {
            throw std::runtime_error("CutManager::addGraphicalCut() - null cut provided!");
        }
        TCutG* owned = dynamic_cast<TCutG*>(cut->Clone());
        owned->SetName(name.c_str());
        graphical_cuts_[name] = GraphicalCut(name, owned, description);
    }
    
    /**
     * @brief Test point against graphical cut
     */
    bool passGraphicalCut(const std::string& name, double x, double y) const {
        auto it = graphical_cuts_.find(name);
        if (it == graphical_cuts_.end()) {
            throw std::runtime_error("CutManager::passGraphicalCut() - Cut '" + name + "' not defined!");
        }
        return it->second.pass(x, y);
    }
    
    /**
     * @brief Check if graphical cut exists
     */
    bool hasGraphicalCut(const std::string& name) const {
        return graphical_cuts_.find(name) != graphical_cuts_.end();
    }
    
    // ========================================================================
    // Cut Activation/Deactivation
    // ========================================================================
    
    void setRangeCutActive(const std::string& name, bool active) {
        getRangeCut(name).active = active;
    }
    
    void setTriggerCutActive(const std::string& name, bool active) {
        auto it = trigger_cuts_.find(name);
        if (it != trigger_cuts_.end()) {
            it->second.active = active;
        }
    }
    
    void setGraphicalCutActive(const std::string& name, bool active) {
        auto it = graphical_cuts_.find(name);
        if (it != graphical_cuts_.end()) {
            it->second.active = active;
        }
    }
    
    void setAllCutsActive(bool active) {
        for (auto& p : range_cuts_) p.second.active = active;
        for (auto& p : trigger_cuts_) p.second.active = active;
        for (auto& p : graphical_cuts_) p.second.active = active;
    }
    
    // ========================================================================
    // Statistics
    // ========================================================================
    
    /**
     * @brief Reset all cut statistics
     */
    void resetStatistics() {
        for (auto& p : range_cuts_) p.second.reset();
        for (auto& p : trigger_cuts_) p.second.reset();
        for (auto& p : graphical_cuts_) p.second.reset();
    }
    
    /**
     * @brief Print cut flow summary
     */
    void printCutFlow(std::ostream& os = std::cout) const {
        os << "\n";
        os << "╔════════════════════════════════════════════════════════════════╗\n";
        os << "║                       CUT FLOW SUMMARY                         ║\n";
        os << "╠════════════════════════════════════════════════════════════════╣\n";
        os << "║ Cut Name                   │ Tested   │ Passed   │ Efficiency  ║\n";
        os << "╠────────────────────────────┼──────────┼──────────┼─────────────╣\n";
        
        // Range cuts
        for (const auto& p : range_cuts_) {
            const auto& c = p.second;
            os << "║ " << std::left << std::setw(26) << c.name 
               << " │ " << std::right << std::setw(8) << c.tested
               << " │ " << std::setw(8) << c.passed
               << " │ " << std::setw(9) << std::fixed << std::setprecision(2) 
               << (c.efficiency() * 100) << "%  ║\n";
        }
        
        // Trigger cuts
        for (const auto& p : trigger_cuts_) {
            const auto& c = p.second;
            os << "║ " << std::left << std::setw(26) << c.name 
               << " │ " << std::right << std::setw(8) << c.tested
               << " │ " << std::setw(8) << c.passed
               << " │ " << std::setw(9) << std::fixed << std::setprecision(2) 
               << (c.efficiency() * 100) << "%  ║\n";
        }
        
        // Graphical cuts
        for (const auto& p : graphical_cuts_) {
            const auto& c = p.second;
            os << "║ " << std::left << std::setw(26) << c.name 
               << " │ " << std::right << std::setw(8) << c.tested
               << " │ " << std::setw(8) << c.passed
               << " │ " << std::setw(9) << std::fixed << std::setprecision(2) 
               << (c.efficiency() * 100) << "%  ║\n";
        }
        
        os << "╚════════════════════════════════════════════════════════════════╝\n";
    }
    
    /**
     * @brief Print defined cuts
     */
    void printDefinedCuts(std::ostream& os = std::cout) const {
        os << "\nDefined Cuts:\n";
        os << "─────────────────────────────────────────\n";
        
        if (!range_cuts_.empty()) {
            os << "Range Cuts:\n";
            for (const auto& p : range_cuts_) {
                os << "  " << p.first << ": [" << p.second.min << ", " << p.second.max << "]";
                if (!p.second.active) os << " (DISABLED)";
                os << "\n";
            }
        }
        
        if (!trigger_cuts_.empty()) {
            os << "Trigger Cuts:\n";
            for (const auto& p : trigger_cuts_) {
                os << "  " << p.first << ": mask=0x" << std::hex << p.second.mask << std::dec;
                if (p.second.require_all) os << " (require ALL)";
                if (!p.second.active) os << " (DISABLED)";
                os << "\n";
            }
        }
        
        if (!graphical_cuts_.empty()) {
            os << "Graphical Cuts:\n";
            for (const auto& p : graphical_cuts_) {
                os << "  " << p.first;
                if (!p.second.active) os << " (DISABLED)";
                os << "\n";
            }
        }
    }
    
    // ========================================================================
    // Query Methods
    // ========================================================================
    
    bool hasRangeCut(const std::string& name) const {
        return range_cuts_.find(name) != range_cuts_.end();
    }
    
    bool hasTriggerCut(const std::string& name) const {
        return trigger_cuts_.find(name) != trigger_cuts_.end();
    }
    
    size_t rangeCutCount() const { return range_cuts_.size(); }
    size_t triggerCutCount() const { return trigger_cuts_.size(); }
    size_t graphicalCutCount() const { return graphical_cuts_.size(); }
    
    std::vector<std::string> listRangeCuts() const {
        std::vector<std::string> names;
        for (const auto& p : range_cuts_) names.push_back(p.first);
        return names;
    }

private:
    std::map<std::string, RangeCut> range_cuts_;
    std::map<std::string, TriggerCut> trigger_cuts_;
    std::map<std::string, GraphicalCut> graphical_cuts_;
};

// ============================================================================
// Convenience Macros for Cut Flow
// ============================================================================

// Helper macro for sequential cuts with early return
#define CUT_FLOW_START(mgr) bool _cutflow_passed = true
#define CUT_FLOW_CHECK(mgr, cutname, value) \
    if (_cutflow_passed && !mgr.passRangeCut(cutname, value)) _cutflow_passed = false
#define CUT_FLOW_CHECK_2D(mgr, cutname, x, y) \
    if (_cutflow_passed && !mgr.passGraphicalCut(cutname, x, y)) _cutflow_passed = false
#define CUT_FLOW_PASSED _cutflow_passed

#endif // CUT_MANAGER_H
