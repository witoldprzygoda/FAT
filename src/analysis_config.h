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
     * @brief Get input source path (can be .root file or .list file)
     * @return Path to input source
     */
    std::string getInputSource() const {
        return config_["input"]["source"].asString();
    }
    
    /**
     * @brief Check if input is a single ROOT file
     * @return true if source ends with .root
     */
    bool isInputRootFile() const {
        std::string source = getInputSource();
        return source.size() >= 5 && 
               source.substr(source.size() - 5) == ".root";
    }
    
    /**
     * @brief Check if input is a file list
     * @return true if source ends with .list
     */
    bool isInputFileList() const {
        std::string source = getInputSource();
        return source.size() >= 5 && 
               source.substr(source.size() - 5) == ".list";
    }
    
    /**
     * @brief Get input file list path (for backward compatibility)
     * @return Path if source is .list file, empty otherwise
     */
    std::string getInputFileList() const {
        return isInputFileList() ? getInputSource() : "";
    }
    
    /**
     * @brief Get input files as vector
     * @return Vector with single .root file or empty if using .list
     */
    std::vector<std::string> getInputFiles() const {
        std::vector<std::string> files;
        if (isInputRootFile()) {
            files.push_back(getInputSource());
        }
        return files;
    }
    
    std::string getInputTreeName() const {
        return config_["input"]["tree_name"].asString("PPip_ID");
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
        os << "╔════════════════════════════════════════════════════════════════╗\n";
        os << "║                   ANALYSIS CONFIGURATION                       ║\n";
        os << "╠════════════════════════════════════════════════════════════════╣\n";
        os << "║ Config file: " << std::left << std::setw(50) << config_file_ << "║\n";
        os << "║                                                                ║\n";
        os << "║ Input:                                                         ║\n";
        std::string source_info = getInputSource();
        if (isInputRootFile()) source_info += " (ROOT file)";
        else if (isInputFileList()) source_info += " (file list)";
        os << "║   Source: " << std::left << std::setw(53) << source_info << "║\n";
        os << "║   Tree: " << std::left << std::setw(55) << getInputTreeName() << "║\n";
        os << "║   Start event: " << std::left << std::setw(48) << getStartEvent() << "║\n";
        os << "║   Max events: " << std::left << std::setw(49) << getMaxEvents() << "║\n";
        os << "║                                                                ║\n";
        os << "║ Output:                                                        ║\n";
        os << "║   File: " << std::left << std::setw(55) << getOutputFilename() << "║\n";
        os << "║                                                                ║\n";
        os << "║ Beam:                                                          ║\n";
        std::ostringstream ke_str;
        ke_str << getBeamKineticEnergy() << " MeV";
        os << "║   Kinetic energy: " << std::left << std::setw(45) << ke_str.str() << "║\n";
        os << "╚════════════════════════════════════════════════════════════════╝\n";
    }

private:
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
