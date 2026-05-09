/**
 * @file ntuple_reader.h
 * @brief Input NTuple reader with reflection-based variable access
 *
 * Provides flexible reading of ROOT TTrees/TNtuples with:
 * - Lazy branch binding (bind on first access)
 * - Named variable access via operator[]
 * - Variable aliasing for channel-transparent analysis
 * - Support for TChain (multiple files)
 * - Automatic type handling for Float_t branches
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef NTUPLE_READER_H
#define NTUPLE_READER_H

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TLeaf.h>
#include <string>
#include <map>
#include <vector>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <fstream>

// ============================================================================
// NTupleReader: Reflection-based input ntuple reader
// ============================================================================
/**
 * @class NTupleReader
 * @brief Reads ROOT TTrees with named variable access
 *
 * Usage Example:
 * @code
 *   NTupleReader reader;
 *   reader.open("data.root", "PPip_ID");
 *   // Or use chain:
 *   reader.openChain({"file1.root", "file2.root"}, "PPip_ID");
 *
 *   for (Long64_t i = 0; i < reader.entries(); ++i) {
 *       reader.getEntry(i);
 *       double p = reader["p_p"];      // Lazy binding
 *       double theta = reader["p_theta"];
 *       // ...
 *   }
 * @endcode
 */
class NTupleReader {
public:
    // ========================================================================
    // Constructors / Destructor
    // ========================================================================
    
    NTupleReader() = default;
    
    ~NTupleReader() {
        // Smart pointers handle cleanup
    }
    
    // Disable copy
    NTupleReader(const NTupleReader&) = delete;
    NTupleReader& operator=(const NTupleReader&) = delete;
    
    // ========================================================================
    // File Opening Methods
    // ========================================================================
    
    /**
     * @brief Open single ROOT file and tree
     * @param filename Path to ROOT file
     * @param treename Name of TTree/TNtuple to read
     */
    void open(const std::string& filename, const std::string& treename) {
        file_ = std::make_unique<TFile>(filename.c_str(), "READ");
        if (!file_ || file_->IsZombie()) {
            throw std::runtime_error("NTupleReader::open() - Cannot open file: " + filename);
        }
        
        tree_ = dynamic_cast<TTree*>(file_->Get(treename.c_str()));
        if (!tree_) {
            throw std::runtime_error("NTupleReader::open() - Tree '" + treename + 
                                   "' not found in " + filename);
        }
        
        treename_ = treename;
        is_chain_ = false;
        // Disable all branches by default; bindBranch() will re-enable each one
        // on first access. Avoids decompressing branches the analysis never reads.
        tree_->SetBranchStatus("*", 0);
        std::cout << "NTupleReader: Opened '" << treename << "' from " << filename
                  << " (" << tree_->GetEntries() << " entries)\n";
    }
    
    /**
     * @brief Open chain of ROOT files
     * @param filenames Vector of file paths
     * @param treename Name of TTree/TNtuple to read
     */
    void openChain(const std::vector<std::string>& filenames, const std::string& treename) {
        chain_ = std::make_unique<TChain>(treename.c_str());
        
        for (const auto& fname : filenames) {
            int added = chain_->Add(fname.c_str());
            if (added == 0) {
                std::cerr << "Warning: No entries added from " << fname << "\n";
            }
        }
        
        if (chain_->GetEntries() == 0) {
            throw std::runtime_error("NTupleReader::openChain() - Chain is empty!");
        }
        
        tree_ = chain_.get();
        treename_ = treename;
        is_chain_ = true;
        // Disable all branches by default; bindBranch() will re-enable each one
        // on first access. SetBranchStatus on a TChain is sticky across files.
        tree_->SetBranchStatus("*", 0);
        std::cout << "NTupleReader: Opened chain '" << treename << "' with "
                  << filenames.size() << " files (" << tree_->GetEntries() << " entries)\n";
    }
    
    /**
     * @brief Open files from a list file
     * @param listfile Path to file containing list of ROOT files
     * @param treename Name of TTree/TNtuple to read
     * 
     * Supports multiple formats:
     * - Plain file paths (one per line)
     * - ROOT code format: chain->Add("/path/to/file.root");
     * - Comments starting with # or //
     */
    void openFromList(const std::string& listfile, const std::string& treename) {
        std::ifstream ifs(listfile);
        if (!ifs) {
            throw std::runtime_error("NTupleReader::openFromList() - Cannot open list file: " + listfile);
        }
        
        std::vector<std::string> files;
        std::string line;
        while (std::getline(ifs, line)) {
            // Skip empty lines
            if (line.empty()) continue;
            
            // Trim leading whitespace
            size_t start = line.find_first_not_of(" \t");
            if (start == std::string::npos) continue;
            
            // Skip comments
            if (line[start] == '#' || line.compare(start, 2, "//") == 0) continue;
            
            std::string filepath;
            
            // Check for chain->Add("...") format
            size_t quote_start = line.find('"');
            if (quote_start != std::string::npos) {
                size_t quote_end = line.find('"', quote_start + 1);
                if (quote_end != std::string::npos) {
                    filepath = line.substr(quote_start + 1, quote_end - quote_start - 1);
                }
            } else {
                // Plain file path format
                size_t end = line.find_last_not_of(" \t;");
                if (end != std::string::npos) {
                    filepath = line.substr(start, end - start + 1);
                }
            }
            
            if (!filepath.empty()) {
                files.push_back(filepath);
            }
        }
        
        if (files.empty()) {
            throw std::runtime_error("NTupleReader::openFromList() - No files in list: " + listfile);
        }
        
        std::cout << "NTupleReader: Found " << files.size() << " files in " << listfile << "\n";
        openChain(files, treename);
    }
    
    // ========================================================================
    // Entry Access
    // ========================================================================
    
    /**
     * @brief Get total number of entries
     */
    Long64_t entries() const {
        if (!tree_) {
            throw std::runtime_error("NTupleReader::entries() - No tree loaded!");
        }
        return tree_->GetEntries();
    }
    
    /**
     * @brief Load specific entry
     * @param entry Entry number to load
     * @return Bytes read (0 if error)
     */
    Int_t getEntry(Long64_t entry) {
        if (!tree_) {
            throw std::runtime_error("NTupleReader::getEntry() - No tree loaded!");
        }
        current_entry_ = entry;
        return tree_->GetEntry(entry);
    }
    
    /**
     * @brief Get current entry number
     */
    Long64_t currentEntry() const {
        return current_entry_;
    }
    
    // ========================================================================
    // Variable Access (Reflection)
    // ========================================================================
    
    /**
     * @brief Access variable by name (lazy binding)
     * @param varname Name of branch/leaf
     * @return Reference to Float_t value
     *
     * On first access, automatically binds the branch.
     * Subsequent accesses use cached address.
     */
    Float_t& operator[](const std::string& varname) {
        auto it = branch_values_.find(varname);
        if (it != branch_values_.end()) {
            return it->second;
        }

        // Resolve prefix mapping: ep_xxx -> prefix1_xxx, em_xxx -> prefix2_xxx
        std::string branch_name = resolveVariable(varname);

        // Lazy binding - bind branch on first access
        // Store under logical name so subsequent lookups are fast
        return bindBranch(varname, branch_name);
    }
    
    /**
     * @brief Const access to variable
     */
    const Float_t& operator[](const std::string& varname) const {
        auto it = branch_values_.find(varname);
        if (it == branch_values_.end()) {
            throw std::runtime_error("NTupleReader::operator[] const - Variable '" + 
                                   varname + "' not bound (use non-const access first)");
        }
        return it->second;
    }
    
    /**
     * @brief Check if variable exists in tree (resolves aliases)
     */
    bool hasVariable(const std::string& varname) const {
        if (!tree_) return false;
        std::string actual = resolveVariable(varname);
        return tree_->GetLeaf(actual.c_str()) != nullptr;
    }
    
    /**
     * @brief Get list of all branch names
     */
    std::vector<std::string> listVariables() const {
        std::vector<std::string> names;
        if (!tree_) return names;
        
        TObjArray* leaves = tree_->GetListOfLeaves();
        for (int i = 0; i < leaves->GetEntries(); ++i) {
            TLeaf* leaf = dynamic_cast<TLeaf*>(leaves->At(i));
            if (leaf) {
                names.push_back(leaf->GetName());
            }
        }
        return names;
    }
    
    /**
     * @brief Pre-bind multiple variables (for performance)
     * @param varnames List of variable names to bind
     */
    void bindVariables(const std::vector<std::string>& varnames) {
        for (const auto& var : varnames) {
            (*this)[var];  // Use operator[] to resolve aliases
        }
    }

    // ========================================================================
    // Variable Prefix Mapping
    // ========================================================================

    /**
     * @brief Set lepton prefix replacements for channel-transparent analysis
     * @param prefix1 Replacement for canonical "ep_" prefix (e.g., "ep1" for EpEp)
     * @param prefix2 Replacement for canonical "em_" prefix (e.g., "ep2" for EpEp)
     *
     * All variables starting with "ep_" will have the prefix replaced with
     * prefix1 + "_", and variables starting with "em_" with prefix2 + "_".
     * Variables with other prefixes (neutr_, fwdet_, isBest, etc.) pass through
     * unchanged.
     *
     * Example:
     * @code
     *   reader.setLeptonPrefixes("ep1", "ep2");  // for EpEp channel
     *   // Now reader["ep_p"] reads from branch "ep1_p"
     *   // and reader["em_p"] reads from branch "ep2_p"
     *   // while reader["neutr_mult"] still reads "neutr_mult"
     * @endcode
     */
    void setLeptonPrefixes(const std::string& prefix1, const std::string& prefix2) {
        lepton_prefix_1_ = prefix1;
        lepton_prefix_2_ = prefix2;
        has_prefix_mapping_ = true;
    }

    /**
     * @brief Check if prefix mapping is active
     */
    bool hasPrefixMapping() const { return has_prefix_mapping_; }
    
    // ========================================================================
    // Accessors
    // ========================================================================
    
    TTree* getTree() { return tree_; }
    const TTree* getTree() const { return tree_; }
    const std::string& getTreeName() const { return treename_; }
    bool isChain() const { return is_chain_; }

    // ------------------------------------------------------------------------
    // Per-file accessors (relevant when isChain() == true): used by the
    // PT3 trigger-bias correction to compute one weight per input file
    // (one HADES run = one ROOT file in the .list).
    //   getCurrentTreeNumber() returns the index of the currently-loaded
    //                          tree in the chain, or 0 for a single TTree.
    //   getNTrees()            returns the number of files in the chain,
    //                          or 1 for a single TTree.
    // ------------------------------------------------------------------------
    int getCurrentTreeNumber() const {
        return is_chain_ ? chain_->GetTreeNumber() : 0;
    }
    int getNTrees() const {
        return is_chain_ ? chain_->GetNtrees() : 1;
    }
    std::string getTreeFilePath(int idx) const {
        if (!is_chain_) {
            return (idx == 0 && file_) ? file_->GetName() : std::string{};
        }
        auto* files = chain_->GetListOfFiles();
        if (!files || idx < 0 || idx >= files->GetEntries()) return {};
        auto* el = files->At(idx);
        return el ? el->GetTitle() : std::string{};
    }

    // ------------------------------------------------------------------------
    // Per-entry / per-file accessors used by trigger calibration: derive the
    // current event's index inside its own tree (chain-aware), and the total
    // entry count of any tree in the chain. Both relative to the input
    // physical file the chain element points at, NOT to the merged chain.
    // ------------------------------------------------------------------------
    Long64_t getLocalEntryInTree() const {
        if (!is_chain_) return current_entry_;
        if (!chain_)    return -1;
        const Long64_t* off = chain_->GetTreeOffset();
        if (!off)       return -1;
        return current_entry_ - off[chain_->GetTreeNumber()];
    }
    Long64_t getTreeNEvents(int idx) const {
        if (!is_chain_) return tree_ ? tree_->GetEntries() : 0;
        if (!chain_ || idx < 0 || idx >= chain_->GetNtrees()) return 0;
        const Long64_t* off = chain_->GetTreeOffset();
        if (!off)       return 0;
        const Long64_t  total = chain_->GetEntries();
        return (idx + 1 < chain_->GetNtrees())
                   ? (off[idx + 1] - off[idx])
                   : (total - off[idx]);
    }
    
    /**
     * @brief Get number of bound variables
     */
    size_t boundVariableCount() const {
        return branch_values_.size();
    }
    
    /**
     * @brief Print summary of reader state
     */
    void printSummary(std::ostream& os = std::cout) const {
        os << "NTupleReader Summary:\n";
        os << "  Tree: " << treename_ << "\n";
        os << "  Type: " << (is_chain_ ? "TChain" : "TTree") << "\n";
        os << "  Entries: " << (tree_ ? tree_->GetEntries() : 0) << "\n";
        os << "  Bound variables: " << branch_values_.size() << "\n";
        if (has_prefix_mapping_) {
            os << "  Prefix mapping: ep_ -> " << lepton_prefix_1_
               << "_, em_ -> " << lepton_prefix_2_ << "_\n";
        }
        if (!branch_values_.empty()) {
            os << "  Variables:\n";
            for (const auto& pair : branch_values_) {
                os << "    - " << pair.first << " = " << pair.second << "\n";
            }
        }
    }

private:
    // ========================================================================
    // Private Methods
    // ========================================================================
    
    /**
     * @brief Bind branch to internal storage
     * @param logical_name Name used for lookup in branch_values_ (what the code uses)
     * @param branch_name Actual branch name in the ROOT tree (may differ if aliased)
     */
    Float_t& bindBranch(const std::string& logical_name, const std::string& branch_name) {
        if (!tree_) {
            throw std::runtime_error("NTupleReader::bindBranch() - No tree loaded!");
        }

        // Check if branch exists (using physical branch name)
        TBranch* branch = tree_->GetBranch(branch_name.c_str());
        if (!branch) {
            // Try leaf (for TNtuple with combined branches)
            TLeaf* leaf = tree_->GetLeaf(branch_name.c_str());
            if (!leaf) {
                std::string msg = "NTupleReader::bindBranch() - Variable '" + logical_name + "'";
                if (logical_name != branch_name) {
                    msg += " (aliased to '" + branch_name + "')";
                }
                msg += " not found in tree '" + treename_ + "'";
                throw std::runtime_error(msg);
            }
        }

        // Create storage under logical name, bind to physical branch.
        // Order matters: re-activate the branch BEFORE SetBranchAddress, otherwise
        // ROOT keeps it skipped during GetEntry and reads return zero.
        branch_values_[logical_name] = 0.0f;
        tree_->SetBranchStatus(branch_name.c_str(), 1);
        tree_->SetBranchAddress(branch_name.c_str(), &branch_values_[logical_name]);

        // Re-read current entry to get value
        if (current_entry_ >= 0) {
            tree_->GetEntry(current_entry_);
        }

        return branch_values_[logical_name];
    }

    /// @brief Convenience overload (no alias)
    Float_t& bindBranch(const std::string& varname) {
        return bindBranch(varname, varname);
    }

    /**
     * @brief Resolve variable name through prefix mapping
     *
     * If prefix mapping is active:
     *   "ep_xxx" -> "prefix1_xxx"  (e.g., "ep_p" -> "ep1_p")
     *   "em_xxx" -> "prefix2_xxx"  (e.g., "em_p" -> "ep2_p")
     * Other variables pass through unchanged (neutr_, fwdet_, isBest, etc.)
     */
    std::string resolveVariable(const std::string& varname) const {
        if (!has_prefix_mapping_) return varname;

        if (varname.compare(0, 3, "ep_") == 0) {
            return lepton_prefix_1_ + varname.substr(2);  // "ep_xxx" -> "prefix1_xxx"
        }
        if (varname.compare(0, 3, "em_") == 0) {
            return lepton_prefix_2_ + varname.substr(2);  // "em_xxx" -> "prefix2_xxx"
        }
        return varname;  // non-lepton variables pass through
    }
    
    // ========================================================================
    // Data Members
    // ========================================================================
    
    std::unique_ptr<TFile> file_;
    std::unique_ptr<TChain> chain_;
    TTree* tree_ = nullptr;  // Points to either file's tree or chain
    
    std::string treename_;
    bool is_chain_ = false;
    Long64_t current_entry_ = -1;
    
    // Storage for branch values (reflection map)
    std::map<std::string, Float_t> branch_values_;

    // Lepton prefix mapping: ep_ -> prefix1_, em_ -> prefix2_
    std::string lepton_prefix_1_ = "ep";   // replacement for "ep" prefix
    std::string lepton_prefix_2_ = "em";   // replacement for "em" prefix
    bool has_prefix_mapping_ = false;
};

#endif // NTUPLE_READER_H
