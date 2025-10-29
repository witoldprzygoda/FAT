# HNtuple Improvements: Enhanced Error Handling and Diagnostics

## Overview

The HNtuple class provides a convenient lazy-construction pattern for ROOT TNtuples, where the structure is determined automatically from the variables you set before the first `fill()` call. This document describes the enhanced error handling and diagnostic features added to make the class more user-friendly and robust.

---

## Key Concept: Lazy Construction with Freeze

### How It Works

```cpp
HNtuple ntuple("data", "My Data");
ntuple.setFile(outFile);

// Define variables BEFORE first fill() - ntuple not yet created
ntuple["energy"] = 100.0;
ntuple["momentum"] = 50.0;
ntuple["theta"] = 45.0;

// First fill() - Creates TNtuple with these 3 variables and FREEZES structure
ntuple.fill();  // ← Structure is now FROZEN!

// Subsequent fills - can modify existing variables but cannot add new ones
ntuple["energy"] = 200.0;  // ✓ OK - modifying existing variable
ntuple["phi"] = 30.0;      // ✗ ERROR - cannot add new variable after freeze!
ntuple.fill();
```

**Key Rule**: All variables must be defined BEFORE the first `fill()` call!

---

## New Features

### 1. Improved Error Messages

#### Before (old error message):
```
A variable: "phi" tried to be assigned after HNtuple was booked,
that is, e.g after the first "fill()" call.
```

#### After (new error message):
```
╔════════════════════════════════════════════════════════════════╗
║  HNtuple ERROR: Cannot add new variable after freeze          ║
╠════════════════════════════════════════════════════════════════╣
║ Attempted to add: "phi"
║ NTuple name:      "data"
║ Fill count:       5 (frozen after fill #1)
║
║ The NTuple structure is FROZEN after the first fill() call.
║ All variables must be defined BEFORE the first fill().
║
║ Current NTuple structure (3 variables):
║ ┌────────────────────────────────────────────────────────────┐
║ │ [0] energy
║ │ [1] momentum
║ │ [2] theta
║ └────────────────────────────────────────────────────────────┘
║
║ SOLUTION:
║   Add 'r["phi"] = value;' BEFORE the first fill() call,
║   or check for typos in the variable name.
╚════════════════════════════════════════════════════════════════╝
```

**Benefits**:
- ✅ Shows exactly which variable you tried to add
- ✅ Lists ALL variables currently in the ntuple
- ✅ Shows fill count (helps debugging)
- ✅ Clear visual formatting with box drawing
- ✅ Provides solution to fix the error
- ✅ Helps catch typos (e.g., "enrgy" vs "energy")

### 2. First Fill() Notification

When the ntuple is frozen (on first `fill()` call), you now get a helpful message:

```
╔════════════════════════════════════════════════════════════════╗
║  HNtuple FROZEN: Structure locked after first fill()          ║
╠════════════════════════════════════════════════════════════════╣
║ NTuple name: "data"
║ Variables:   3
║ ┌────────────────────────────────────────────────────────────┐
║ │ [0] energy
║ │ [1] momentum
║ │ [2] theta
║ └────────────────────────────────────────────────────────────┘
║
║ This structure is now FROZEN. No new variables can be added.
╚════════════════════════════════════════════════════════════════╝
```

**Benefits**:
- ✅ Confirms the structure was frozen
- ✅ Shows exactly what variables are in your ntuple
- ✅ Reminder that no new variables can be added

### 3. Validation: Empty NTuple Check

If you try to `fill()` without setting any variables, you get a clear error:

```cpp
HNtuple ntuple("empty", "Empty NTuple");
ntuple.setFile(outFile);
ntuple.fill();  // ERROR!
```

**Error message**:
```
HNtuple ERROR: Attempting to fill() without setting any variables!
NTuple "empty" has no variables defined.
Use: r["variable_name"] = value; before calling fill().
```

### 4. New Query Methods

#### `isFrozen()` - Check if structure is frozen
```cpp
if (ntuple.isFrozen()) {
    std::cout << "Structure is frozen - cannot add variables\n";
} else {
    std::cout << "Can still add variables\n";
}
```

#### `getNVariables()` - Get number of variables
```cpp
std::cout << "NTuple has " << ntuple.getNVariables() << " variables\n";
```

#### `getVariableNames()` - Get list of variable names
```cpp
auto names = ntuple.getVariableNames();  // Returns std::vector<std::string>
for (const auto& name : names) {
    std::cout << "  - " << name << "\n";
}
```

#### `hasVariable(key)` - Check if variable exists
```cpp
if (ntuple.hasVariable("energy")) {
    std::cout << "NTuple has 'energy' variable\n";
}
```

### 5. New Diagnostic Methods

#### `printStructure()` - Print full structure to console
```cpp
ntuple.printStructure();  // Prints to std::cout

// Or to a file
std::ofstream logfile("ntuple_structure.log");
ntuple.printStructure(logfile);
```

**Output**:
```
╔════════════════════════════════════════════════════════════════╗
║  HNtuple Structure                                             ║
╠════════════════════════════════════════════════════════════════╣
║ Name:        "data"
║ Title:       "My Data"
║ Status:      FROZEN
║ Fill count:  42
║ Variables:   3
║ ┌────────────────────────────────────────────────────────────┐
║ │ [0] energy = 100.5
║ │ [1] momentum = 50.2
║ │ [2] theta = 45.0
║ └────────────────────────────────────────────────────────────┘
╚════════════════════════════════════════════════════════════════╝
```

#### `getStructureString()` - Get structure as string
```cpp
std::string structure = ntuple.getStructureString();
// Store in database, log file, or send to monitoring system
```

### 6. Fill Count Tracking

The ntuple now tracks how many times `fill()` has been called:

```cpp
// Internal member (accessible via printStructure())
Long64_t fillCount;  // Increments on each fill()
```

This helps with debugging and understanding when the freeze happened.

---

## Common Use Cases

### Use Case 1: Checking Structure Before Processing

```cpp
void analyzeData(HNtuple& ntuple) {
    // Check if ntuple has the variables we need
    if (!ntuple.hasVariable("energy") || !ntuple.hasVariable("momentum")) {
        std::cerr << "ERROR: Required variables missing!\n";
        ntuple.printStructure(std::cerr);
        return;
    }

    // Proceed with analysis
    // ...
}
```

### Use Case 2: Debugging Variable Name Typos

```cpp
// You write:
ntuple["enrgy"] = 100.0;  // Typo!

// You get a detailed error showing all valid names:
// Current NTuple structure:
//   [0] energy    ← Notice the correct spelling!
//   [1] momentum
//   [2] theta
```

The error message makes it immediately obvious that you typed "enrgy" instead of "energy".

### Use Case 3: Dynamic Variable Discovery

```cpp
void processNTuple(HNtuple& ntuple) {
    std::cout << "Processing ntuple with variables:\n";
    auto variables = ntuple.getVariableNames();
    for (const auto& var : variables) {
        std::cout << "  - " << var << "\n";
    }
}
```

### Use Case 4: Conditional Variable Addition

```cpp
HNtuple ntuple("data", "Data");
ntuple.setFile(outFile);

// Add basic variables
ntuple["energy"] = 100.0;
ntuple["momentum"] = 50.0;

// Conditionally add simulation variables
if (isMC) {
    ntuple["sim_energy"] = 98.5;
    ntuple["sim_momentum"] = 49.8;
}

// First fill - structure is now frozen with 2 or 4 variables
ntuple.fill();

// Later code can check:
if (ntuple.hasVariable("sim_energy")) {
    // This is MC data
}
```

### Use Case 5: Logging NTuple Structure

```cpp
void logAnalysisSetup(HNtuple& ntuple, const std::string& logfile) {
    std::ofstream log(logfile, std::ios::app);
    log << "Analysis started: " << std::time(nullptr) << "\n";
    log << "NTuple structure:\n";
    ntuple.printStructure(log);
    log << "\n";
}
```

---

## Migration from Old Code

### Old Code (no error handling):
```cpp
HNtuple ntuple("data", "Data");
ntuple.setFile(outFile);

ntuple["energy"] = 100.0;
ntuple["momentum"] = 50.0;
ntuple.fill();

// Later, you forget structure is frozen:
ntuple["phi"] = 30.0;  // Cryptic error message
ntuple.fill();
```

### New Code (with better errors):
```cpp
HNtuple ntuple("data", "Data");
ntuple.setFile(outFile);

ntuple["energy"] = 100.0;
ntuple["momentum"] = 50.0;
ntuple.fill();  // ← Prints structure and confirms freeze

// Later, when you try to add a new variable:
ntuple["phi"] = 30.0;  // ← Clear, detailed error message showing:
                       //   - What you tried to add ("phi")
                       //   - What variables exist (energy, momentum)
                       //   - How to fix it
```

**No code changes needed** - your existing code works the same, but errors are much more informative!

---

## Technical Details

### Header Changes (`hntuple.h`)

**Added members**:
```cpp
Long64_t fillCount {0};  // Track number of fills
```

**Added methods**:
```cpp
Bool_t isFrozen() const;
Int_t getNVariables() const;
std::vector<std::string> getVariableNames() const;
Bool_t hasVariable(const std::string& key) const;
void printStructure(std::ostream& os = std::cout) const;
std::string getStructureString() const;
```

### Implementation Changes (`hntuple.cc`)

1. **`operator[]` (line 117-146)**: Enhanced error message with structure display
2. **`const operator[]` (line 161-173)**: Improved error for reading non-existent variables
3. **`fill()` (line 177-259)**:
   - Added validation (empty ntuple check)
   - Added freeze notification
   - Added fillCount increment
4. **New methods** (line 261-322): Implementations of query/diagnostic methods

---

## Performance Impact

**Zero performance impact** in normal operation:
- Error handling only executes when errors occur
- `fillCount` increment is O(1)
- Query methods are only called when explicitly requested
- No overhead during regular `fill()` operations

---

## Testing

Run the comprehensive test suite:

```bash
cd /home/user/FAT/examples
make test-hntuple
```

This runs 5 examples demonstrating:
1. ✅ Correct usage (lazy construction)
2. ✅ Error when adding variable after freeze
3. ✅ New query/diagnostic methods
4. ✅ Error when filling without variables
5. ✅ Typo detection in variable names

---

## Summary

### Problems Solved

1. **❌ Old**: Cryptic error messages
   **✅ New**: Detailed, actionable error messages

2. **❌ Old**: No way to query ntuple state
   **✅ New**: Multiple query methods (isFrozen, hasVariable, etc.)

3. **❌ Old**: No structure visibility
   **✅ New**: printStructure() shows full structure with current values

4. **❌ Old**: Typos hard to debug
   **✅ New**: Error shows all valid variables, making typos obvious

5. **❌ Old**: No validation
   **✅ New**: Prevents filling empty ntuples

### Key Benefits

- ✅ **Better debugging**: Clear error messages with context
- ✅ **Safer code**: Validation catches mistakes early
- ✅ **More maintainable**: Query methods for introspection
- ✅ **User-friendly**: Visual formatting with box drawing
- ✅ **Production-ready**: Logging support via printStructure()
- ✅ **Backward compatible**: Existing code works unchanged

---

## Example Error Message Comparison

### Before
```
terminate called after throwing an instance of 'std::invalid_argument'
  what():  A variable: "phi" tried to be assigned after HNtuple was booked
```

### After
```
╔════════════════════════════════════════════════════════════════╗
║  HNtuple ERROR: Cannot add new variable after freeze          ║
╠════════════════════════════════════════════════════════════════╣
║ Attempted to add: "phi"
║ NTuple name:      "data"
║ Fill count:       5 (frozen after fill #1)
║
║ Current NTuple structure (3 variables):
║ ┌────────────────────────────────────────────────────────────┐
║ │ [0] energy
║ │ [1] momentum
║ │ [2] theta
║ └────────────────────────────────────────────────────────────┘
║
║ SOLUTION:
║   Add 'r["phi"] = value;' BEFORE the first fill() call
╚════════════════════════════════════════════════════════════════╝
```

**The difference is night and day!** 🌟

---

## Conclusion

The improved HNtuple error handling transforms debugging from frustrating guesswork into straightforward problem-solving. The enhanced error messages immediately show you:
- What went wrong
- What the current state is
- How to fix it

This makes the lazy-construction pattern much more user-friendly and production-ready!
