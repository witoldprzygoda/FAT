/**
 * @file ntuple_reader.h
 * @brief Input NTuple reader with reflection-based variable access
 *
 * Provides flexible reading of ROOT TTrees/TNtuples with:
 * - Lazy branch binding (bind on first access)
 * - Named variable access via operator[]
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
        
        // Lazy binding - bind branch on first access
        return bindBranch(varname);
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
     * @brief Check if variable exists in tree
     */
    bool hasVariable(const std::string& varname) const {
        if (!tree_) return false;
        return tree_->GetLeaf(varname.c_str()) != nullptr;
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
            bindBranch(var);
        }
    }
    
    // ========================================================================
    // Accessors
    // ========================================================================
    
    TTree* getTree() { return tree_; }
    const TTree* getTree() const { return tree_; }
    const std::string& getTreeName() const { return treename_; }
    bool isChain() const { return is_chain_; }
    
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
     */
    Float_t& bindBranch(const std::string& varname) {
        if (!tree_) {
            throw std::runtime_error("NTupleReader::bindBranch() - No tree loaded!");
        }
        
        // Check if branch exists
        TBranch* branch = tree_->GetBranch(varname.c_str());
        if (!branch) {
            // Try leaf (for TNtuple with combined branches)
            TLeaf* leaf = tree_->GetLeaf(varname.c_str());
            if (!leaf) {
                throw std::runtime_error("NTupleReader::bindBranch() - Variable '" + 
                                       varname + "' not found in tree '" + treename_ + "'");
            }
        }
        
        // Create storage and bind
        branch_values_[varname] = 0.0f;
        tree_->SetBranchAddress(varname.c_str(), &branch_values_[varname]);
        
        // Re-read current entry to get value
        if (current_entry_ >= 0) {
            tree_->GetEntry(current_entry_);
        }
        
        return branch_values_[varname];
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
};

#endif // NTUPLE_READER_H
