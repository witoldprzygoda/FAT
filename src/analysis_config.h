/**
 * @file analysis_config.h
 * @brief JSON-based analysis configuration system
 *
 * Provides external configuration for:
 * - Input files (lists, chains)
 * - NTuple names
 * - Output file settings
 * - Cut definitions
 * - Beam parameters
 * - Analysis flags
 *
 * Uses a simple built-in JSON parser (no external dependencies).
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef ANALYSIS_CONFIG_H
#define ANALYSIS_CONFIG_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cctype>

// Forward declare KinematicType (defined in pparticle.h)
enum class KinematicType;

// Forward declare ConsoleBox (defined in console_box.h)
class ConsoleBox;

// ============================================================================
// Simple JSON Value (lightweight implementation)
// ============================================================================

class JsonValue {
public:
    enum Type { NONE, BOOL, NUMBER, STRING, ARRAY, OBJECT };
    
    JsonValue() : type_(NONE) {}
    JsonValue(bool b) : type_(BOOL), bool_val_(b) {}
    JsonValue(double n) : type_(NUMBER), num_val_(n) {}
    JsonValue(int n) : type_(NUMBER), num_val_(n) {}
    JsonValue(const std::string& s) : type_(STRING), str_val_(s) {}
    JsonValue(const char* s) : type_(STRING), str_val_(s) {}
    
    Type type() const { return type_; }
    bool isNull() const { return type_ == NONE; }
    bool isBool() const { return type_ == BOOL; }
    bool isNumber() const { return type_ == NUMBER; }
    bool isString() const { return type_ == STRING; }
    bool isArray() const { return type_ == ARRAY; }
    bool isObject() const { return type_ == OBJECT; }
    
    // Getters with defaults
    bool asBool(bool def = false) const {
        return type_ == BOOL ? bool_val_ : def;
    }
    
    double asDouble(double def = 0.0) const {
        return type_ == NUMBER ? num_val_ : def;
    }
    
    int asInt(int def = 0) const {
        return type_ == NUMBER ? static_cast<int>(num_val_) : def;
    }
    
    std::string asString(const std::string& def = "") const {
        return type_ == STRING ? str_val_ : def;
    }
    
    // Array access
    const std::vector<JsonValue>& asArray() const {
        static std::vector<JsonValue> empty;
        return type_ == ARRAY ? arr_val_ : empty;
    }
    
    size_t size() const {
        if (type_ == ARRAY) return arr_val_.size();
        if (type_ == OBJECT) return obj_val_.size();
        return 0;
    }
    
    const JsonValue& operator[](size_t index) const {
        static JsonValue null_val;
        if (type_ != ARRAY || index >= arr_val_.size()) return null_val;
        return arr_val_[index];
    }
    
    // Object access
    const JsonValue& operator[](const std::string& key) const {
        static JsonValue null_val;
        if (type_ != OBJECT) return null_val;
        auto it = obj_val_.find(key);
        return it != obj_val_.end() ? it->second : null_val;
    }
    
    bool has(const std::string& key) const {
        if (type_ != OBJECT) return false;
        return obj_val_.find(key) != obj_val_.end();
    }
    
    std::vector<std::string> keys() const {
        std::vector<std::string> k;
        if (type_ == OBJECT) {
            for (const auto& p : obj_val_) k.push_back(p.first);
        }
        return k;
    }
    
    // Builders
    static JsonValue array() {
        JsonValue v;
        v.type_ = ARRAY;
        return v;
    }
    
    static JsonValue object() {
        JsonValue v;
        v.type_ = OBJECT;
        return v;
    }
    
    void push_back(const JsonValue& v) {
        if (type_ == NONE) type_ = ARRAY;
        if (type_ == ARRAY) arr_val_.push_back(v);
    }
    
    void set(const std::string& key, const JsonValue& v) {
        if (type_ == NONE) type_ = OBJECT;
        if (type_ == OBJECT) obj_val_[key] = v;
    }
    
private:
    Type type_;
    bool bool_val_ = false;
    double num_val_ = 0.0;
    std::string str_val_;
    std::vector<JsonValue> arr_val_;
    std::map<std::string, JsonValue> obj_val_;
};

// ============================================================================
// Simple JSON Parser
// ============================================================================

class JsonParser {
public:
    static JsonValue parse(const std::string& json) {
        size_t pos = 0;
        return parseValue(json, pos);
    }
    
    static JsonValue parseFile(const std::string& filename) {
        std::ifstream ifs(filename);
        if (!ifs) {
            throw std::runtime_error("JsonParser::parseFile() - Cannot open: " + filename);
        }
        
        std::stringstream ss;
        ss << ifs.rdbuf();
        return parse(ss.str());
    }

private:
    static void skipWhitespace(const std::string& s, size_t& pos) {
        while (pos < s.size() && std::isspace(s[pos])) ++pos;
        // Skip comments (// style)
        if (pos + 1 < s.size() && s[pos] == '/' && s[pos+1] == '/') {
            while (pos < s.size() && s[pos] != '\n') ++pos;
            skipWhitespace(s, pos);
        }
    }
    
    static JsonValue parseValue(const std::string& s, size_t& pos) {
        skipWhitespace(s, pos);
        if (pos >= s.size()) return JsonValue();
        
        char c = s[pos];
        
        if (c == '{') return parseObject(s, pos);
        if (c == '[') return parseArray(s, pos);
        if (c == '"') return parseString(s, pos);
        if (c == 't' || c == 'f') return parseBool(s, pos);
        if (c == 'n') return parseNull(s, pos);
        if (std::isdigit(c) || c == '-' || c == '+') return parseNumber(s, pos);
        
        throw std::runtime_error("JsonParser: Unexpected character at position " + 
                               std::to_string(pos));
    }
    
    static JsonValue parseObject(const std::string& s, size_t& pos) {
        JsonValue obj = JsonValue::object();
        ++pos;  // skip '{'
        
        skipWhitespace(s, pos);
        if (s[pos] == '}') {
            ++pos;
            return obj;
        }
        
        while (true) {
            skipWhitespace(s, pos);
            
            // Parse key
            if (s[pos] != '"') {
                throw std::runtime_error("JsonParser: Expected string key at position " + 
                                       std::to_string(pos));
            }
            std::string key = parseString(s, pos).asString();
            
            skipWhitespace(s, pos);
            if (s[pos] != ':') {
                throw std::runtime_error("JsonParser: Expected ':' at position " + 
                                       std::to_string(pos));
            }
            ++pos;
            
            // Parse value
            JsonValue val = parseValue(s, pos);
            obj.set(key, val);
            
            skipWhitespace(s, pos);
            if (s[pos] == '}') {
                ++pos;
                return obj;
            }
            if (s[pos] != ',') {
                throw std::runtime_error("JsonParser: Expected ',' or '}' at position " + 
                                       std::to_string(pos));
            }
            ++pos;
        }
    }
    
    static JsonValue parseArray(const std::string& s, size_t& pos) {
        JsonValue arr = JsonValue::array();
        ++pos;  // skip '['
        
        skipWhitespace(s, pos);
        if (s[pos] == ']') {
            ++pos;
            return arr;
        }
        
        while (true) {
            JsonValue val = parseValue(s, pos);
            arr.push_back(val);
            
            skipWhitespace(s, pos);
            if (s[pos] == ']') {
                ++pos;
                return arr;
            }
            if (s[pos] != ',') {
                throw std::runtime_error("JsonParser: Expected ',' or ']' at position " + 
                                       std::to_string(pos));
            }
            ++pos;
        }
    }
    
    static JsonValue parseString(const std::string& s, size_t& pos) {
        ++pos;  // skip opening '"'
        std::string result;
        
        while (pos < s.size() && s[pos] != '"') {
            if (s[pos] == '\\' && pos + 1 < s.size()) {
                ++pos;
                switch (s[pos]) {
                    case 'n': result += '\n'; break;
                    case 't': result += '\t'; break;
                    case 'r': result += '\r'; break;
                    case '"': result += '"'; break;
                    case '\\': result += '\\'; break;
                    default: result += s[pos];
                }
            } else {
                result += s[pos];
            }
            ++pos;
        }
        
        if (pos >= s.size()) {
            throw std::runtime_error("JsonParser: Unterminated string");
        }
        ++pos;  // skip closing '"'
        return JsonValue(result);
    }
    
    static JsonValue parseNumber(const std::string& s, size_t& pos) {
        size_t start = pos;
        if (s[pos] == '-' || s[pos] == '+') ++pos;
        
        while (pos < s.size() && (std::isdigit(s[pos]) || s[pos] == '.' || 
               s[pos] == 'e' || s[pos] == 'E' || s[pos] == '+' || s[pos] == '-')) {
            ++pos;
        }
        
        std::string num_str = s.substr(start, pos - start);
        return JsonValue(std::stod(num_str));
    }
    
    static JsonValue parseBool(const std::string& s, size_t& pos) {
        if (s.compare(pos, 4, "true") == 0) {
            pos += 4;
            return JsonValue(true);
        }
        if (s.compare(pos, 5, "false") == 0) {
            pos += 5;
            return JsonValue(false);
        }
        throw std::runtime_error("JsonParser: Invalid boolean at position " + 
                               std::to_string(pos));
    }
    
    static JsonValue parseNull(const std::string& s, size_t& pos) {
        if (s.compare(pos, 4, "null") == 0) {
            pos += 4;
            return JsonValue();
        }
        throw std::runtime_error("JsonParser: Invalid null at position " + 
                               std::to_string(pos));
    }
};

// ============================================================================
// AnalysisConfig: Main configuration class
// ============================================================================

/**
 * @class AnalysisConfig
 * @brief Loads and provides access to analysis configuration
 *
 * Configuration file format (JSON):
 * @code
 * {
 *   "input": {
 *     "file_list": "h68_10.list",
 *     "tree_name": "PPip_ID",
 *     "max_events": -1
 *   },
 *   "output": {
 *     "filename": "output.root"
 *   },
 *   "beam": {
 *     "kinetic_energy": 1580.0
 *   },
 *   "cuts": {
 *     "neutron_mass": {"min": 0.899, "max": 0.986},
 *     "deltaPP_mass": {"min": 0.8, "max": 1.8}
 *   },
 *   "triggers": {
 *     "physics": {"mask": 4, "require_all": false}
 *   },
 *   "graphical_cuts": {
 *     "proton_pid": {"file": "cuts/proton.root", "name": "proton_cut"}
 *   }
 * }
 * @endcode
 */
class AnalysisConfig {
public:
    // ========================================================================
    // Loading
    // ========================================================================
    
    /**
     * @brief Load configuration from JSON file
     */
    void load(const std::string& filename) {
        config_ = JsonParser::parseFile(filename);
        config_file_ = filename;
        std::cout << "AnalysisConfig: Loaded configuration from " << filename << "\n";
        
        // Validate required fields
        validate();
    }
    
    /**
     * @brief Load configuration from JSON string
     */
    void loadFromString(const std::string& json) {
        config_ = JsonParser::parse(json);
        config_file_ = "<string>";
        validate();
    }
    
    // ========================================================================
    // Input Configuration
    // ========================================================================
    
    /**
     * @brief Get input source path (can be .root file, .list file, or comma-separated .root files)
     * @return Path to input source (first file if multiple)
     *
     * Supported formats:
     * - Single ROOT file: "source": "file.root"
     * - File list: "source": "files.list"
     * - Comma-separated: "source": "file1.root, file2.root, file3.root"
     * - JSON array: "source": ["file1.root", "file2.root"]
     */
    std::string getInputSource() const {
        // Check if source is a JSON array
        const JsonValue& source_val = config_["input"]["source"];
        if (source_val.isArray() && source_val.size() > 0) {
            return source_val[0].asString();
        }
        return source_val.asString();
    }

    /**
     * @brief Check if input is ROOT file(s) - single, comma-separated, or array
     * @return true if any source contains .root
     */
    bool isInputRootFile() const {
        const JsonValue& source_val = config_["input"]["source"];

        // JSON array of ROOT files
        if (source_val.isArray()) {
            if (source_val.size() > 0) {
                std::string first = source_val[0].asString();
                return first.size() >= 5 && first.substr(first.size() - 5) == ".root";
            }
            return false;
        }

        // String - could be single file, list, or comma-separated
        std::string source = source_val.asString();

        // Check for comma-separated (contains comma and .root)
        if (source.find(',') != std::string::npos && source.find(".root") != std::string::npos) {
            return true;
        }

        // Single file ending with .root
        return source.size() >= 5 && source.substr(source.size() - 5) == ".root";
    }

    /**
     * @brief Check if input is a file list (.list file)
     * @return true if source ends with .list
     */
    bool isInputFileList() const {
        const JsonValue& source_val = config_["input"]["source"];
        if (source_val.isArray()) return false;

        std::string source = source_val.asString();
        return source.size() >= 5 &&
               source.substr(source.size() - 5) == ".list";
    }

    /**
     * @brief Check if input is multiple ROOT files (comma-separated or array)
     * @return true if multiple ROOT files specified
     */
    bool isInputMultipleRootFiles() const {
        const JsonValue& source_val = config_["input"]["source"];

        // JSON array with multiple files
        if (source_val.isArray() && source_val.size() > 1) {
            return true;
        }

        // Comma-separated string
        std::string source = source_val.asString();
        return source.find(',') != std::string::npos && source.find(".root") != std::string::npos;
    }

    /**
     * @brief Get input file list path (for backward compatibility)
     * @return Path if source is .list file, empty otherwise
     */
    std::string getInputFileList() const {
        return isInputFileList() ? getInputSource() : "";
    }

    /**
     * @brief Get all input ROOT files as vector
     * @return Vector of ROOT file paths (handles single, comma-separated, and array formats)
     *
     * Examples:
     * - "file.root" → ["file.root"]
     * - "a.root, b.root" → ["a.root", "b.root"]
     * - ["a.root", "b.root"] → ["a.root", "b.root"]
     */
    std::vector<std::string> getInputFiles() const {
        std::vector<std::string> files;
        const JsonValue& source_val = config_["input"]["source"];

        // JSON array format
        if (source_val.isArray()) {
            for (size_t i = 0; i < source_val.size(); ++i) {
                std::string file = source_val[i].asString();
                // Trim whitespace
                file.erase(0, file.find_first_not_of(" \t"));
                file.erase(file.find_last_not_of(" \t") + 1);
                if (!file.empty()) {
                    files.push_back(file);
                }
            }
            return files;
        }

        // String format
        std::string source = source_val.asString();

        // Check for comma-separated files
        if (source.find(',') != std::string::npos && source.find(".root") != std::string::npos) {
            std::stringstream ss(source);
            std::string file;
            while (std::getline(ss, file, ',')) {
                // Trim whitespace
                file.erase(0, file.find_first_not_of(" \t"));
                file.erase(file.find_last_not_of(" \t") + 1);
                if (!file.empty()) {
                    files.push_back(file);
                }
            }
            return files;
        }

        // Single ROOT file
        if (isInputRootFile()) {
            files.push_back(source);
        }

        return files;
    }
    
    std::string getInputTreeName() const {
        return config_["input"]["tree_name"].asString("PPip_ID");
    }

    /**
     * @brief Get channel name (e.g., "EpEm", "EpEp", "EmEm")
     * @return Channel name from config, empty string if not set
     */
    std::string getChannelName() const {
        return config_["input"]["channel"].asString("");
    }

    /**
     * @brief Get lepton prefix replacements for channel-transparent analysis
     * @return Pair of {prefix1, prefix2} replacing canonical "ep"/"em" prefixes
     *
     * Used to run the same analysis code on different trees (EpEm, EpEp, EmEm).
     * In config: "lepton_prefixes": ["ep1", "ep2"]  (for EpEp channel)
     * If not set, returns empty pair (no prefix mapping needed).
     */
    std::pair<std::string, std::string> getLeptonPrefixes() const {
        const JsonValue& arr = config_["input"]["lepton_prefixes"];
        if (arr.isArray() && arr.size() >= 2) {
            return {arr[0].asString(), arr[1].asString()};
        }
        return {"", ""};
    }

    /**
     * @brief Get starting event number
     * @return Starting event index (default 0)
     */
    Long64_t getStartEvent() const {
        return static_cast<Long64_t>(config_["input"]["start_event"].asDouble(0));
    }
    
    Long64_t getMaxEvents() const {
        return static_cast<Long64_t>(config_["input"]["max_events"].asDouble(-1));
    }
    
    // ========================================================================
    // Output Configuration
    // ========================================================================
    
    std::string getOutputFilename() const {
        return config_["output"]["filename"].asString("output.root");
    }
    
    std::string getOutputOption() const {
        return config_["output"]["option"].asString("RECREATE");
    }
    
    /**
     * @brief Get whether to keep intermediate TTree file
     * @return true to keep, false to delete after conversion (default: false)
     */
    bool getKeepIntermediateTree() const {
        return config_["output"]["keep_intermediate_tree"].asBool(false);
    }
    
    /**
     * @brief Get missing value for DynamicHNtuple
     * @return Sentinel value for missing variables (default: -1.0)
     */
    Float_t getMissingValue() const {
        return static_cast<Float_t>(config_["output"]["missing_value"].asDouble(-1.0));
    }
    
    // ========================================================================
    // Beam Configuration
    // ========================================================================
    
    double getBeamKineticEnergy() const {
        return config_["beam"]["kinetic_energy"].asDouble(1580.0);
    }
    
    /**
     * @brief Get luminosity scaling factor for histograms
     * @return Luminosity factor (default: 1.0 = no scaling)
     * 
     * Used to normalize histograms to absolute cross-section units.
     * Example: If luminosity = 2.5e-3, all histogram fills will be scaled
     * by this factor when using weighted fill methods.
     */
    double getLuminosity() const {
        return config_["beam"]["luminosity"].asDouble(1.0);
    }
    
    // ========================================================================
    // Kinematics Configuration
    // ========================================================================
    
    /**
     * @brief Check if RECONSTRUCTED kinematic type is enabled
     * @return true if reconstructed data should be loaded (default: true)
     */
    bool hasReconstructed() const {
        return config_["kinematics"]["reconstructed"].asBool(true);
    }
    
    /**
     * @brief Check if CORRECTED kinematic type is enabled
     * @return true if corrected data should be loaded (default: true)
     */
    bool hasCorrected() const {
        return config_["kinematics"]["corrected"].asBool(true);
    }
    
    /**
     * @brief Check if SIMULATED kinematic type is enabled
     * @return true if simulated data should be loaded (default: false)
     */
    bool hasSimulated() const {
        return config_["kinematics"]["simulated"].asBool(false);
    }
    
    /**
     * @brief Get the kinematic type to use for analysis
     * 
     * Returns the type specified in config.json under kinematics.analysis_type.
     * Valid values: "reconstructed", "corrected", "simulated"
     * Default: RECONSTRUCTED
     * 
     * @return KinematicType enum value
     */
    KinematicType getAnalysisKinematicType() const {
        std::string type_str = config_["kinematics"]["analysis_type"].asString("reconstructed");
        
        // Convert to lowercase for comparison
        std::string lower;
        for (char c : type_str) {
            lower += std::tolower(static_cast<unsigned char>(c));
        }
        
        if (lower == "corrected") {
            return KinematicType::CORRECTED;
        } else if (lower == "simulated") {
            return KinematicType::SIMULATED;
        }
        // Default to RECONSTRUCTED
        return KinematicType::RECONSTRUCTED;
    }
    
    // ========================================================================
    // Forward Tracker Configuration
    // ========================================================================
    
    /**
     * @brief Check if Forward Tracker processing is enabled
     * @return true if FWD objects should be built (default: false)
     */
    bool isFwdEnabled() const {
        return config_["fwdet"].asBool(false);
    }

    // ========================================================================
    // ECAL (Electromagnetic Calorimeter) Configuration
    // ========================================================================

    /**
     * @brief Check if ECAL processing is enabled
     * @return true if ECAL objects should be built (default: false)
     */
    bool isEcalEnabled() const {
        return config_["ecal"].asBool(false);
    }

    // ========================================================================
    // Trigger
    //
    // The full-mode analysis no longer cuts on the trigbit: every event with
    // trigbit ∈ {8192 (PT3), 4096 (PT2)} is processed, and Manager auto-routes
    // every fill to either the histogram (PT3) or its "_pt2" twin (PT2). The
    // bin-by-bin PT2/PT3 ratio is the trigger correction. No JSON knobs for
    // trigger selection or bias correction — the cal-mode is the only path
    // that still treats the trigbit specially (and it accepts everything).
    // ========================================================================

    // ========================================================================
    // Run mode
    // ========================================================================
    //   "full" (default) — full analysis: histograms, dilepton_nt, ecal_nt,
    //                      etc. and the existing Pass 1 trigger-bias counts.
    //   "trigger_calibration" — minimal pass that only fills trigger_cal_nt
    //                      (raw trigbit + oa + isBest + vertex + start
    //                       per event, no cuts at fill time) plus a small
    //                      trigger_cal_files metadata TTree. Pass 1, all
    //                      other ntuples and all histograms are skipped.
    //                      Used by research/calibration/ to derive the
    //                      per-segment PT3-bias weights.
    std::string getMode() const {
        return config_["mode"].asString("full");
    }
    bool isTriggerCalibrationOnly() const {
        return getMode() == "trigger_calibration";
    }

    // ========================================================================
    // Cut Configuration
    // ========================================================================
    
    /**
     * @brief Get all defined range cuts
     * @return Map of cut name -> {min, max}
     */
    std::map<std::string, std::pair<double, double>> getRangeCuts() const {
        std::map<std::string, std::pair<double, double>> cuts;
        const JsonValue& cuts_obj = config_["cuts"];
        
        for (const auto& name : cuts_obj.keys()) {
            const JsonValue& cut = cuts_obj[name];
            double min_val = cut["min"].asDouble(0);
            double max_val = cut["max"].asDouble(0);
            cuts[name] = {min_val, max_val};
        }
        return cuts;
    }
    
    /**
     * @brief Get trigger cut definitions
     */
    struct TriggerDef {
        int mask;
        bool require_all;
    };
    
    std::map<std::string, TriggerDef> getTriggerCuts() const {
        std::map<std::string, TriggerDef> triggers;
        const JsonValue& trig_obj = config_["triggers"];
        
        for (const auto& name : trig_obj.keys()) {
            const JsonValue& t = trig_obj[name];
            TriggerDef def;
            def.mask = t["mask"].asInt(0);
            def.require_all = t["require_all"].asBool(false);
            triggers[name] = def;
        }
        return triggers;
    }
    
    /**
     * @brief Get graphical cut definitions
     */
    struct GraphicalCutDef {
        std::string file;
        std::string name;  // Name in ROOT file (optional, defaults to key)
    };
    
    std::map<std::string, GraphicalCutDef> getGraphicalCuts() const {
        std::map<std::string, GraphicalCutDef> gcuts;
        const JsonValue& gc_obj = config_["graphical_cuts"];
        
        for (const auto& key : gc_obj.keys()) {
            const JsonValue& g = gc_obj[key];
            GraphicalCutDef def;
            def.file = g["file"].asString();
            def.name = g["name"].asString(key);  // Default to key name
            gcuts[key] = def;
        }
        return gcuts;
    }
    
    // ========================================================================
    // Variable Configuration (optional - which variables to read)
    // ========================================================================
    
    std::vector<std::string> getRequiredVariables() const {
        std::vector<std::string> vars;
        const JsonValue& arr = config_["variables"];
        for (size_t i = 0; i < arr.size(); ++i) {
            vars.push_back(arr[i].asString());
        }
        return vars;
    }
    
    // ========================================================================
    // Custom Parameters
    // ========================================================================
    
    /**
     * @brief Get custom parameter by path (e.g., "analysis.weight_mode")
     */
    JsonValue get(const std::string& path) const {
        std::vector<std::string> parts;
        std::string part;
        for (char c : path) {
            if (c == '.') {
                if (!part.empty()) {
                    parts.push_back(part);
                    part.clear();
                }
            } else {
                part += c;
            }
        }
        if (!part.empty()) parts.push_back(part);
        
        const JsonValue* current = &config_;
        for (const auto& p : parts) {
            current = &((*current)[p]);
            if (current->isNull()) break;
        }
        return *current;
    }
    
    double getDouble(const std::string& path, double def = 0.0) const {
        return get(path).asDouble(def);
    }
    
    int getInt(const std::string& path, int def = 0) const {
        return get(path).asInt(def);
    }
    
    std::string getString(const std::string& path, const std::string& def = "") const {
        return get(path).asString(def);
    }
    
    bool getBool(const std::string& path, bool def = false) const {
        return get(path).asBool(def);
    }
    
    // ========================================================================
    // Diagnostics
    // ========================================================================
    
    void print(std::ostream& os = std::cout) const {
        os << "\n";
        
        // Title
        os << "╔════════════════════════════════════════════════════════════════╗\n";
        os << "║                   ANALYSIS CONFIGURATION                       ║\n";
        os << "╠════════════════════════════════════════════════════════════════╣\n";
        
        // Config file
        os << "║ Config file: " << std::left << std::setw(50) << config_file_ << "║\n";
        os << "║                                                                ║\n";
        
        // Input section
        os << "║ Input:                                                         ║\n";
        if (isInputMultipleRootFiles()) {
            std::vector<std::string> files = getInputFiles();
            std::string source_info = std::to_string(files.size()) + " ROOT files (chain)";
            printConfigLine(os, "Source", source_info);
        } else {
            std::string source_info = getInputSource();
            if (isInputRootFile()) source_info += " (ROOT file)";
            else if (isInputFileList()) source_info += " (file list)";
            printConfigLine(os, "Source", source_info);
        }
        printConfigLine(os, "Tree", getInputTreeName());
        printConfigLine(os, "Start event", std::to_string(getStartEvent()));
        printConfigLine(os, "Max events", std::to_string(getMaxEvents()));
        os << "║                                                                ║\n";
        
        // Output section
        os << "║ Output:                                                        ║\n";
        printConfigLine(os, "File", getOutputFilename());
        os << "║                                                                ║\n";
        
        // Beam section
        os << "║ Beam:                                                          ║\n";
        std::ostringstream ke_str;
        ke_str << getBeamKineticEnergy() << " MeV";
        printConfigLine(os, "Kinetic energy", ke_str.str());
        os << "║                                                                ║\n";

        // Trigger handling is now implicit: PT3 → original histograms,
        // PT2 → "_pt2" twins, both written to the output ROOT file.
        os << "║ Trigger:  auto-route per event (PT3 → H, PT2 → H_pt2)          ║\n";
        os << "╚════════════════════════════════════════════════════════════════╝\n";
    }

private:
    static constexpr int BOX_INNER_WIDTH = 64;  // Inner width of config box
    
    void printConfigLine(std::ostream& os, const std::string& label, 
                        const std::string& value) const {
        std::string line = "  " + label + ": " + value;
        int padding = BOX_INNER_WIDTH - static_cast<int>(line.length());
        if (padding < 0) padding = 0;
        os << "║" << line;
        for (int i = 0; i < padding; ++i) os << ' ';
        os << "║\n";
    }
    
    void validate() const {
        // Check for required sections
        if (!config_.has("input")) {
            std::cerr << "Warning: Config missing 'input' section\n";
        }
        if (!config_.has("output")) {
            std::cerr << "Warning: Config missing 'output' section\n";
        }
    }
    
    JsonValue config_;
    std::string config_file_;
};

#endif // ANALYSIS_CONFIG_H
