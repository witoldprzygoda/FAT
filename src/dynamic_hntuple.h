/**
 * @file dynamic_hntuple.h
 * @brief Dynamic HNtuple with unlimited variable discovery
 *
 * This class allows adding variables at ANY time during processing.
 * Internally uses TTree for flexible storage, then converts to
 * flat TNtuple at finalization for easy plotting.
 *
 * Key features:
 * - Add variables at any time via operator[]
 * - Intermediate TTree storage (handles dynamic schema)
 * - Final conversion to flat TNtuple (alphabetically ordered)
 * - Missing values filled with configurable sentinel (default: -1)
 * - Progress indicator during conversion
 *
 * Usage:
 * @code
 *   DynamicHNtuple nt("name", "title", output_file);
 *   
 *   // In event loop - add ANY variable at ANY time
 *   nt["event"] = event_num;
 *   nt["mass"] = calculated_mass;
 *   if (has_hit_3) {
 *       nt["hit_3"] = hit_value;  // OK even if first time!
 *   }
 *   nt.fill();
 *   
 *   // At the end
 *   nt.finalize();  // Converts TTree -> TNtuple
 * @endcode
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef DYNAMIC_HNTUPLE_H
#define DYNAMIC_HNTUPLE_H

#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <cstdio>
#include <chrono>
#include <fstream>

// ============================================================================
// DynamicHNtuple: Unlimited variable discovery with TTree→TNtuple conversion
// ============================================================================

class DynamicHNtuple {
public:
    // ========================================================================
    // Constructor / Destructor
    // ========================================================================
    
    /**
     * @brief Construct DynamicHNtuple
     * @param name Name of the ntuple
     * @param title Title/description
     * @param output_file Pointer to final output ROOT file
     * @param missing_value Value to use for missing variables (default: -1)
     * @param keep_intermediate Keep intermediate TTree file (default: false)
     */
    DynamicHNtuple(const std::string& name, const std::string& title,
                   TFile* output_file,
                   Float_t missing_value = -1.0f,
                   bool keep_intermediate = false)
        : name_(name), title_(title), output_file_(output_file),
          missing_value_(missing_value), keep_intermediate_(keep_intermediate)
    {
        if (!output_file_) {
            throw std::runtime_error("DynamicHNtuple: output_file cannot be null!");
        }
        
        // Create intermediate file with unique name per ntuple
        // Format: output_name_tree.root (e.g., output_ppip_nt_particles_tree.root)
        std::string out_path = output_file_->GetName();
        size_t dot_pos = out_path.rfind('.');
        if (dot_pos != std::string::npos) {
            intermediate_filename_ = out_path.substr(0, dot_pos) + "_" + name_ + "_tree.root";
        } else {
            intermediate_filename_ = out_path + "_" + name_ + "_tree.root";
        }
        
        // Open intermediate file
        intermediate_file_ = std::make_unique<TFile>(intermediate_filename_.c_str(), "RECREATE");
        if (!intermediate_file_ || intermediate_file_->IsZombie()) {
            throw std::runtime_error("DynamicHNtuple: Cannot create intermediate file: " + intermediate_filename_);
        }
        
        // Create TTree
        tree_ = new TTree((name_ + "_tree").c_str(), title_.c_str());
        tree_->SetDirectory(intermediate_file_.get());
        
        std::cout << "DynamicHNtuple: Created '" << name_ << "' with intermediate storage: " 
                  << intermediate_filename_ << "\n";
    }
    
    ~DynamicHNtuple() {
        // Clean up branch value storage
        for (auto& pair : branch_values_) {
            delete pair.second;
        }
    }
    
    // Disable copy
    DynamicHNtuple(const DynamicHNtuple&) = delete;
    DynamicHNtuple& operator=(const DynamicHNtuple&) = delete;
    
    // ========================================================================
    // Variable Access - Add ANY variable at ANY time
    // ========================================================================
    
    /**
     * @brief Access/create variable by name
     * 
     * If variable doesn't exist, it's created automatically.
     * Safe to call at any point during processing.
     */
    Float_t& operator[](const std::string& key) {
        if (finalized_) {
            throw std::runtime_error("DynamicHNtuple: Cannot add variables after finalize()!");
        }
        
        auto it = branch_values_.find(key);
        if (it != branch_values_.end()) {
            return *(it->second);
        }
        
        // New variable - create branch dynamically
        Float_t* value_ptr = new Float_t(missing_value_);
        branch_values_[key] = value_ptr;
        discovered_vars_.insert(key);
        
        // Create branch in TTree
        tree_->Branch(key.c_str(), value_ptr, (key + "/F").c_str());
        
        return *value_ptr;
    }
    
    /**
     * @brief Const access (for reading)
     */
    const Float_t& operator[](const std::string& key) const {
        auto it = branch_values_.find(key);
        if (it != branch_values_.end()) {
            return *(it->second);
        }
        throw std::runtime_error("DynamicHNtuple: Variable '" + key + "' not found!");
    }
    
    // ========================================================================
    // Fill - Store current event
    // ========================================================================
    
    /**
     * @brief Fill current event to TTree
     * 
     * After fill(), all values are reset to missing_value for next event.
     */
    void fill() {
        if (finalized_) {
            throw std::runtime_error("DynamicHNtuple: Cannot fill() after finalize()!");
        }
        
        tree_->Fill();
        fill_count_++;
        
        // Reset all values to missing for next event
        for (auto& pair : branch_values_) {
            *(pair.second) = missing_value_;
        }
    }
    
    // ========================================================================
    // Finalize - Convert TTree to TNtuple
    // ========================================================================
    
    /**
     * @brief Finalize: Convert intermediate TTree to flat TNtuple
     * 
     * This must be called before closing the output file.
     * Converts all data from TTree to TNtuple format with:
     * - All discovered variables (alphabetically ordered)
     * - Missing values filled with configured sentinel
     * - Progress indicator during conversion
     */
    void finalize() {
        if (finalized_) {
            std::cerr << "Warning: DynamicHNtuple::finalize() called multiple times\n";
            return;
        }
        
        // Handle case with no variables or no entries
        if (discovered_vars_.empty() || fill_count_ == 0) {
            std::cout << "DynamicHNtuple '" << name_ << "': ";
            if (discovered_vars_.empty()) {
                std::cout << "No variables defined.\n";
            } else {
                std::cout << "No entries filled (" << discovered_vars_.size() << " variables defined).\n";
            }
            
            // Clean up intermediate file
            cleanupIntermediateFile();
            
            // Create empty TNtuple if variables were defined
            if (!discovered_vars_.empty()) {
                std::vector<std::string> sorted_vars(discovered_vars_.begin(), discovered_vars_.end());
                std::string varlist;
                for (size_t i = 0; i < sorted_vars.size(); ++i) {
                    varlist += sorted_vars[i];
                    if (i < sorted_vars.size() - 1) varlist += ":";
                }
                output_file_->cd();
                auto ntuple = std::make_unique<TNtuple>(name_.c_str(), title_.c_str(), varlist.c_str());
                ntuple->Write();
                std::cout << "✓ Created empty TNtuple '" << name_ << "' with " << sorted_vars.size() << " variables\n";
            }
            
            finalized_ = true;
            return;
        }
        
        std::cout << "\n";
        std::cout << "╔════════════════════════════════════════════════════════════════╗\n";
        std::cout << "║         Converting TTree → TNtuple                             ║\n";
        std::cout << "╠════════════════════════════════════════════════════════════════╣\n";
        std::cout << "║ NTuple:     " << std::left << std::setw(51) << name_ << " ║\n";
        std::cout << "║ Variables:  " << std::left << std::setw(51) << discovered_vars_.size() << " ║\n";
        std::cout << "║ Entries:    " << std::left << std::setw(51) << fill_count_ << " ║\n";
        std::cout << "╚════════════════════════════════════════════════════════════════╝\n";
        
        // Build variable list (alphabetical order from std::set)
        std::vector<std::string> sorted_vars(discovered_vars_.begin(), discovered_vars_.end());
        
        std::string varlist;
        for (size_t i = 0; i < sorted_vars.size(); ++i) {
            varlist += sorted_vars[i];
            if (i < sorted_vars.size() - 1) varlist += ":";
        }
        
        // Print variable list
        std::cout << "\nVariables (alphabetical order):\n";
        for (size_t i = 0; i < sorted_vars.size(); ++i) {
            std::cout << "  [" << std::setw(2) << i << "] " << sorted_vars[i] << "\n";
        }
        std::cout << "\n";
        
        // Write and close intermediate tree
        intermediate_file_->cd();
        tree_->Write();
        
        // Reopen for reading
        intermediate_file_->Close();
        intermediate_file_ = std::make_unique<TFile>(intermediate_filename_.c_str(), "READ");
        TTree* read_tree = dynamic_cast<TTree*>(intermediate_file_->Get((name_ + "_tree").c_str()));
        
        if (!read_tree) {
            throw std::runtime_error("DynamicHNtuple: Failed to reopen intermediate TTree!");
        }
        
        // Set up branch addresses for reading
        std::map<std::string, Float_t> read_values;
        for (const auto& var : sorted_vars) {
            read_values[var] = missing_value_;
            TBranch* branch = read_tree->GetBranch(var.c_str());
            if (branch) {
                read_tree->SetBranchAddress(var.c_str(), &read_values[var]);
            }
        }
        
        // Create final TNtuple in output file
        output_file_->cd();
        auto ntuple = std::make_unique<TNtuple>(name_.c_str(), title_.c_str(), varlist.c_str());
        
        // Prepare float array for filling
        std::vector<Float_t> values(sorted_vars.size());
        
        // Convert with progress indicator
        Long64_t total = read_tree->GetEntries();
        Long64_t last_percent = -1;
        
        std::cout << "Converting: ";
        std::cout.flush();
        
        auto start_time = std::chrono::steady_clock::now();
        
        for (Long64_t i = 0; i < total; ++i) {
            // Reset to missing value before reading
            for (auto& pair : read_values) {
                pair.second = missing_value_;
            }
            
            read_tree->GetEntry(i);
            
            // Copy values in order
            for (size_t j = 0; j < sorted_vars.size(); ++j) {
                values[j] = read_values[sorted_vars[j]];
            }
            
            ntuple->Fill(values.data());
            
            // Progress indicator
            Long64_t percent = (i + 1) * 100 / total;
            if (percent != last_percent) {
                last_percent = percent;
                
                // Calculate ETA
                auto now = std::chrono::steady_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
                
                int bar_width = 40;
                int filled = static_cast<int>(percent * bar_width / 100);
                
                std::cout << "\r[";
                for (int b = 0; b < bar_width; ++b) {
                    if (b < filled) std::cout << "█";
                    else std::cout << "░";
                }
                std::cout << "] " << std::setw(3) << percent << "%";
                
                if (percent > 0 && percent < 100) {
                    Long64_t eta = elapsed * (100 - percent) / percent;
                    std::cout << "  ETA: " << std::setw(2) << eta / 60 << ":" 
                              << std::setw(2) << std::setfill('0') << eta % 60 
                              << std::setfill(' ');
                }
                std::cout.flush();
            }
        }
        
        auto end_time = std::chrono::steady_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
        
        std::cout << "\r[";
        for (int b = 0; b < 40; ++b) std::cout << "█";
        std::cout << "] 100%  Done in " << total_time / 60 << ":" 
                  << std::setw(2) << std::setfill('0') << total_time % 60 
                  << std::setfill(' ') << "     \n";
        
        // Write ntuple to output file
        output_file_->cd();
        ntuple->Write();
        
        // Close intermediate file
        intermediate_file_->Close();
        intermediate_file_.reset();
        
        // Clean up intermediate file
        cleanupIntermediateFile();
        
        finalized_ = true;
        
        std::cout << "✓ TNtuple '" << name_ << "' created with " << sorted_vars.size() 
                  << " variables, " << total << " entries\n";
    }
    
    // ========================================================================
    // Query Methods
    // ========================================================================
    
    const std::string& getName() const { return name_; }
    const std::string& getTitle() const { return title_; }
    
    bool isFinalized() const { return finalized_; }
    Long64_t getFillCount() const { return fill_count_; }
    size_t getVariableCount() const { return discovered_vars_.size(); }
    
    std::vector<std::string> getVariableNames() const {
        return std::vector<std::string>(discovered_vars_.begin(), discovered_vars_.end());
    }
    
    bool hasVariable(const std::string& name) const {
        return discovered_vars_.find(name) != discovered_vars_.end();
    }
    
    /**
     * @brief Clean up intermediate file
     */
    void cleanupIntermediateFile() {
        // Close intermediate file if still open
        if (intermediate_file_) {
            intermediate_file_->Close();
            intermediate_file_.reset();
        }
        
        // Delete file unless keeping
        if (!keep_intermediate_) {
            if (std::remove(intermediate_filename_.c_str()) == 0) {
                std::cout << "✓ Removed intermediate file: " << intermediate_filename_ << "\n";
            } else {
                // File might not exist or already deleted - that's OK
                // Only warn if file exists but couldn't be deleted
                std::ifstream test(intermediate_filename_);
                if (test.good()) {
                    std::cerr << "Warning: Could not remove intermediate file: " << intermediate_filename_ << "\n";
                }
            }
        } else {
            std::cout << "✓ Kept intermediate file: " << intermediate_filename_ << "\n";
        }
    }
    
    /**
     * @brief Print current structure
     */
    void printStructure(std::ostream& os = std::cout) const {
        os << "DynamicHNtuple '" << name_ << "':\n";
        os << "  Status: " << (finalized_ ? "FINALIZED" : "COLLECTING") << "\n";
        os << "  Fill count: " << fill_count_ << "\n";
        os << "  Variables (" << discovered_vars_.size() << "):\n";
        int idx = 0;
        for (const auto& var : discovered_vars_) {
            os << "    [" << idx++ << "] " << var << "\n";
        }
    }

private:
    std::string name_;
    std::string title_;
    TFile* output_file_;  // Final output file (not owned)
    
    std::string intermediate_filename_;
    std::unique_ptr<TFile> intermediate_file_;
    TTree* tree_ = nullptr;  // Owned by intermediate_file_
    
    std::map<std::string, Float_t*> branch_values_;  // Current event values
    std::set<std::string> discovered_vars_;          // All discovered variable names (sorted)
    
    Float_t missing_value_ = -1.0f;
    bool keep_intermediate_ = false;
    bool finalized_ = false;
    Long64_t fill_count_ = 0;
};

#endif // DYNAMIC_HNTUPLE_H
